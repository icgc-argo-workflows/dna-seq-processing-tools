#!/usr/bin/env python3

import os
import sys
import json
import re
import glob
from argparse import ArgumentParser
import zipfile
import hashlib
import copy

"""
Major steps:
- convert input Seq to unaligned BAM for each read group
"""


def calculate_size(file_path):
    return os.stat(file_path).st_size

def calculate_md5(file_path):
    md5 = hashlib.md5()
    with open(file_path, 'rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            md5.update(chunk)
    return md5.hexdigest()

def get_files_info(file_to_upload):
    payload_files = {}
    payload_files['name'] = os.path.basename(file_to_upload)
    payload_files['local_path'] = file_to_upload
    payload_files['size'] = calculate_size(file_to_upload)
    payload_files['checksum'] = calculate_md5(file_to_upload)

    return payload_files

def main(args):
    cwd = os.getcwd()
    with zipfile.ZipFile(args.payload_schema_zip) as myzip:
        myzip.extractall(cwd)

    payload_template = glob.glob("argo-metadata-schemas-*/schemas/_example_docs/*.%s.*.json" % args.payload_type)

    output = {}
    if args.payload_type == 'lane_seq_submission':
        with open(args.metadata_lane_seq, 'r') as f:
            metadata = json.load(f)
        if args.input_seq_format == 'FASTQ':
            read_group = metadata.get("read_groups")
            for template in payload_template:
                if not '36.lane_seq_submission.01.ok.json' in template: continue
                with open(template, 'r') as f:
                    payload = json.load(f)

            payload['program'] = metadata.get('program')

            #get inputs of the payload
            for rg in read_group:
                rg_id = rg.get("submitter_id")
                rg_fname = "".join([c if re.match(r"[a-zA-Z0-9\-_]", c) else "_" for c in rg_id])
                if not rg_fname in args.file_to_upload: continue
                payload['inputs']['read_group_submitter_id'] = rg_id
                payload['inputs']['files']['fastq'] = rg.get('files')

        elif args.input_seq_format == 'BAM':
            files = metadata.get("files")
            for template in payload_template:
                if not '35.lane_seq_submission.01.ok.json' in template: continue
                with open(template, 'r') as f:
                    payload = json.load(f)

            payload['program'] = metadata.get('program')

            # get inputs of the payload
            for input_file in files:
                for rg in input_file.get('read_groups'):
                    rg_id = rg.get("submitter_id")
                    rg_fname = "".join([c if re.match(r"[a-zA-Z0-9\-_]", c) else "_" for c in rg_id])
                    if not rg_fname in args.file_to_upload: continue
                    payload['inputs']['read_group_submitter_id'] = rg_id
                    payload['inputs']['files']['bam'] = copy.deepcopy(input_file)
                    payload['inputs']['files']['bam'].pop('read_groups')

        else:
            sys.exit('\n%s: Input files format are not FASTQ or BAM')

        #get files of the payload
        payload['files']['bam_file'].update(get_files_info(args.file_to_upload))

        payload['files']['bam_file'].pop('_final_doc', None)
        payload['files']['bam_file'].pop('_mocked_system_properties', None)

    elif args.payload_type == 'dna_alignment':
        for template in payload_template:
            if not '40.dna_alignment.01.ok.json' in template: continue
            with open(template, 'r') as f:
                payload = json.load(f)

        lane_seq_list = []
        for res_file in args.lane_seq_analysis:
            lane_seq = {}
            with open(res_file, 'r') as f:
                res_json = json.load(f)
            payload['program'] = res_json.get('program')

            lane_seq['lane_seq_submission_id'] = res_json.get('id')
            lane_seq['files'] = {}
            lane_seq['files']['lane_seq'] = res_json['files']['bam_file']

            lane_seq['files']['lane_seq'].update({"bundle_id": res_json.get('id')})

            lane_seq_list.append(lane_seq)

        payload['inputs']['lane_seq'] = lane_seq_list

        #get files of the payload
        payload['files']['aligned_seq'].update(get_files_info(args.file_to_upload))

        #get index files of the payload
        if not os.path.exists(args.file_to_upload + ".bai"): sys.exit('\n%s: Missing BAI index file')
        payload['files']['aligned_seq_index'].update(get_files_info(args.file_to_upload + ".bai"))

        payload['files']['aligned_seq'].pop('_final_doc', None)
        payload['files']['aligned_seq'].pop('_mocked_system_properties', None)
        payload['files']['aligned_seq_index'].pop('_final_doc', None)
        payload['files']['aligned_seq_index'].pop('_mocked_system_properties', None)

    elif args.payload_type == 'sanger_ssm_call':
        pass

    else:
        sys.exit('\n%s: Unknown payload_type')


    payload.pop('_final_doc', None)
    payload.pop('_mocked_system_properties', None)

    # get analysis of the payload
    pass

    payload_fname = ".".join([args.payload_type, os.path.basename(args.file_to_upload), 'json'])
    with open(payload_fname, 'w') as f:
        f.write(json.dumps(payload))

    # write the parameter to stdout
    print(json.dumps(output))

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-sf", "--input_seq_format", dest="input_seq_format",
                        help="Sequence format")
    parser.add_argument("-pt", "--payload_type", dest="payload_type",
                        help="Payload type")
    parser.add_argument("-ps", "--payload_schema_zip", dest="payload_schema_zip", help="released metadata schema zip file")
    parser.add_argument("-mls", "--metadata_lane_seq", dest="metadata_lane_seq",
                        help="json file containing experiment, read_group and file information for sequence preprocessing")
    parser.add_argument("-fu", "--file_to_upload", dest="file_to_upload", help="File to upload to server")
    parser.add_argument("-lp", "--lane_seq_analysis", dest="lane_seq_analysis", help="Analysis of lane seq submission",
                        type=str, nargs='+')
    args = parser.parse_args()

    main(args)