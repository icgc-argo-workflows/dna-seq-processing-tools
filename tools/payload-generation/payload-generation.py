#!/usr/bin/env python3

import os
import sys
import json
import re
from argparse import ArgumentParser
import hashlib
import copy
import requests

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

def run_cmd(cmd):
    print('command: %s' % cmd)
    stdout, stderr, p, success = '', '', None, True
    try:
        p = subprocess.Popen([cmd],
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             shell=True)
        stdout, stderr = p.communicate()
    except Exception as e:
        print('Execution failed: %s' % e, file=sys.stderr)
        success = False

    if p and p.returncode != 0:
        print('Execution failed, none zero code returned.', file=sys.stderr)
        success = False

    print(stdout.decode("utf-8"))
    print(stderr.decode("utf-8"), file=sys.stderr)

    if not success:
        sys.exit(p.returncode if p.returncode else 1)

    return


def main(args):

    if args.bundle_type == 'lane_seq_submission':
        with open(args.metadata_lane_seq, 'r') as f:
            metadata = json.load(f)
        if args.input_seq_format == 'FASTQ':
            read_group = metadata.get("read_groups")

            payload_template_url = "https://raw.githubusercontent.com/icgc-argo/argo-metadata-schemas/%s/schemas/_example_docs/36.lane_seq_submission.01.ok.json" % args.payload_schema_version
            cmd = "curl -sSL -o template --retry 10 %s" % payload_template_url
            run_cmd(cmd)

            with open("template", "r") as f:
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

            payload_template_url = "https://raw.githubusercontent.com/icgc-argo/argo-metadata-schemas/%s/schemas/_example_docs/35.lane_seq_submission.01.ok.json" % args.payload_schema_version
            cmd = "curl -sSL -o template --retry 10 %s" % payload_template_url
            run_cmd(cmd)

            with open("template", "r") as f:
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

    elif args.bundle_type == 'dna_alignment':
        payload_template_url = "https://raw.githubusercontent.com/icgc-argo/argo-metadata-schemas/%s/schemas/_example_docs/40.dna_alignment.01.ok.json" % args.payload_schema_version
        cmd = "curl -sSL -o template --retry 10 %s" % payload_template_url
        run_cmd(cmd)

        with open("template", "r") as f:
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
        if os.path.exists(args.file_to_upload + ".bai"):
            payload['files']['aligned_seq_index'].update(get_files_info(args.file_to_upload + ".bai"))
        elif os.path.exists(args.file_to_upload + ".crai"):
            payload['files']['aligned_seq_index'].update(get_files_info(args.file_to_upload + ".crai"))
        else:
            sys.exit('\n%s: Missing index file')

        payload['files']['aligned_seq'].pop('_final_doc', None)
        payload['files']['aligned_seq'].pop('_mocked_system_properties', None)
        payload['files']['aligned_seq_index'].pop('_final_doc', None)
        payload['files']['aligned_seq_index'].pop('_mocked_system_properties', None)

    elif args.bundle_type == 'sanger_ssm_call':
        pass

    else:
        sys.exit('\n%s: Unknown bundle_type')


    payload.pop('_final_doc', None)
    payload.pop('_mocked_system_properties', None)

    # get analysis of the payload
    pass

    payload_fname = ".".join([args.bundle_type, os.path.basename(args.file_to_upload), 'json'])
    with open(payload_fname, 'w') as f:
        f.write(json.dumps(payload))

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-s", "--input_seq_format", dest="input_seq_format",
                        help="Sequence format")
    parser.add_argument("-t", "--bundle_type", dest="bundle_type",
                        help="Payload type")
    parser.add_argument("-p", "--payload_schema_version", dest="payload_schema_version", help="release version of payload schema")
    parser.add_argument("-m", "--metadata_lane_seq", dest="metadata_lane_seq",
                        help="json file containing experiment, read_group and file information for sequence preprocessing")
    parser.add_argument("-f", "--file_to_upload", dest="file_to_upload", help="File to upload to server")
    parser.add_argument("-a", "--lane_seq_analysis", dest="lane_seq_analysis", help="Analysis of lane seq submission",
                        type=str, nargs='+')
    args = parser.parse_args()

    main(args)
