#!/usr/bin/env python3

import sys
import os
from argparse import ArgumentParser
import json
import subprocess

def bam_check(args):
    with open(args.metadata_json, 'r') as f:
        metadata = json.load(f)

    files = metadata.get('files')

    for _file in files:
        file_with_path = os.path.join(args.seq_files_dir, _file.get('name'))

        # check whether the file exist
        if not os.path.isfile(file_with_path): sys.exit('\n The file: %s does not exist!' % file_with_path)

        rg_metadata = set()
        for rg in _file.get('read_groups'):
            rg_metadata.add(rg.get('submitter_id'))

        # retrieve the @RG from BAM header
        try:
            header = subprocess.check_output(['samtools', 'view', '-H', file_with_path])

        except Exception as e:
            sys.exit('\n%s: Retrieve BAM header failed: %s' % (e, file_with_path))

        # get @RG
        header_array = header.decode('utf-8').rstrip().split('\n')
        rg_bam = set()
        for line in header_array:
            if not line.startswith("@RG"): continue
            rg_array = line.rstrip().split('\t')[1:]
            for element in rg_array:
                if not element.startswith("ID"): continue
                rg_bam.add(element.replace("ID:", ""))

        # compare the RG ids
        if not rg_metadata == rg_bam:
            sys.exit('\nThe read group Ids in metadata do not match with those in BAM file: %s!' % file_with_path)  # die fast

    return True


def fastq_check(args):
    return True


def run_validation(args):
    output_json = {}

    if args.seq_format == "BAM":
        valid = bam_check(args)
    elif args.seq_format == "FASTQ":
        valid = fastq_check(args)
    else:
        sys.exit('\nError: The input files should have format in BAM or FASTQ.')

    if not valid:
        sys.exit('\nError: The input files failed at the validation.')

    output_json['valid'] = 'valid'

    # write the parameter to stdout
    print(json.dumps(output_json))

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-i", "--seq_format", dest="seq_format",
                        help="Sequence format")
    parser.add_argument("-d", "--seq_files_dir", dest="seq_files_dir", help="Directory with seq files to submit and process")
    parser.add_argument("-p", "--metadata_json", dest="metadata_json",
                        help="json file containing experiment, read_group and file information for sequence preprocessing")
    args = parser.parse_args()


    run_validation(args)