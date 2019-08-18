#!/usr/bin/env python3

import os
import sys
import csv
from argparse import ArgumentParser
import subprocess


def run_command(cmd):
    p = subprocess.Popen(cmd,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE,
                         shell=True)
    stdout, stderr = p.communicate()

    return (p, stdout, stderr)


def main(args):
    if args.seq_files:  # if provided as local file
        return  # do nothing
    else:
        if not (args.files_tsv and args.repository and args.token_file):
            sys.exit('Error, when seq-files is missing must provide files-tsv, repository and token-file')

        # read access token and set ENV variable
        with open(args.token_file, 'r') as f:
            os.environ["ACCESSTOKEN"] = f.read().strip()

        files_to_download = set()
        # read TSV to get all unique files for download
        with open(args.files_tsv) as f:
            rd = csv.DictReader(f, delimiter="\t")
            for row in rd:
                files_to_download.add((row['local_path'], row['name']))

        for f in files_to_download:
            file_uri, file_name = f
            if not file_uri.startswith('score://'):
                sys.exit("Only support file_uri with 'score://' scheme, actual file_uri: %s" % file_uri)

            _, _, repo, _, object_id = file_uri.split('/')
            if repo != args.repository:
                sys.exit("File URI repository: %s differs from what's specified: %s" %
                            (repo, args.repository)
                        )

            ret = run_command('score-client --profile %s download --object-id %s --output-dir .' %
                                 (repo, object_id)
                             )

            if ret[0].returncode != 0:
                sys.exit('Download failed on: %s. Error msg: %s' % (file_uri, ret[2]))

            if not os.path.isfile(file_name):
                sys.exit("Object '%s' downloaded but it does not match the expected file name '%s'" % 
                            (file_uri, file_name)
                        )


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-s", "--seq-files", dest="seq_files", nargs="+")
    parser.add_argument("-f", "--files-tsv", dest="files_tsv")
    parser.add_argument("-r", "--repository", dest="repository", choices=['collab','aws'])
    parser.add_argument("-t", "--token-file", dest="token_file")
    args = parser.parse_args()

    main(args)
