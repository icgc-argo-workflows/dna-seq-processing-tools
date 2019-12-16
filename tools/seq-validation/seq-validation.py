#!/usr/bin/env python3

"""
 Copyright (c) 2019, Ontario Institute for Cancer Research (OICR).
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU Affero General Public License as published
 by the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Affero General Public License for more details.
 You should have received a copy of the GNU Affero General Public License
 along with this program. If not, see <https://www.gnu.org/licenses/>.
 Authors:
   Linda Xiang <linda.xiang@oicr.on.ca>
   Junjun Zhang <junjun.zhang@oicr.on.ca>
 """

import sys
import os
from argparse import ArgumentParser
import json
import subprocess
import hashlib


def run_cmd(cmd):
    p = subprocess.Popen(cmd,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE,
                         shell=True)
    stdout, stderr = p.communicate()

    return (p, stdout, stderr)


def calculate_md5(file_path):
    md5 = hashlib.md5()
    with open(file_path, 'rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            md5.update(chunk)
    return md5.hexdigest()


def picard_validate_bam(bam_file):
    # to be implemented
    pass


def get_rg_from_bam_header(bam_file, rg_info):
    cmd = "samtools view -H %s |grep '^@RG'" % bam_file
    p, stdout, stderr = run_cmd(cmd)
    if p.returncode != 0:
        sys.exit("Unable to run 'samtools' to get header, BAM file: %s. Error: %s" % (bam_file, stderr))

    for r in stdout.decode("utf-8").splitlines():
        rg = {}
        for kv in r.split('\t'):
            if ':' not in kv: continue
            tokens = kv.split(':')
            if tokens[0] in rg:
                sys.exit("@RG header error, TAG '%s' occurred more than once in '%s' from BAM file: %s" % \
                    (tokens[0], r, bam_file))

            rg[tokens[0]] = ':'.join(tokens[1:])

        rg_id = rg.get('ID')
        if not rg_id: sys.exit("@RG header error, ID not defined in '%s' from BAM file: %s" % (r, bam_file))

        rg['file'] = bam_file

        if rg_id in rg_info: # same RG ID appeared before
            sys.exit("Same RG ID (%s) exists more than once. First in '%s', second in '%s'" % \
                (rg_id, rg_info[rg_id]['file'], rg['file']))

        rg_info[rg_id] = rg


def rg_info_validation(metadata, rg_info):
    # to be implemented
    print(rg_info)


def file_check(metadata, seq_files):
    file_name_to_path = {}
    for f in seq_files:
        file_name_to_path[os.path.basename(f)] = f

    for f in metadata.get('files'):
        file_name = f.get('name')

        # file existence check: a.1
        if not os.path.isfile(file_name_to_path.get(file_name, "")):
            sys.exit("Specified file does not exist: %s" % file_name)

        # file size check: a.2
        actual_file_size = os.path.getsize(file_name_to_path.get(file_name))
        if actual_file_size != f.get('size'):
            sys.exit("Size of file '%s' does not match what specified in file TSV: %s vs %s" % \
                (file_name, actual_file_size, f.get('size')))

        # file md5sum check: a.3
        actual_file_md5sum = calculate_md5(file_name_to_path.get(file_name))
        if actual_file_md5sum != f.get('md5sum'):
            sys.exit("Md5sum of file '%s' does not match what specified in file TSV: %s vs %s" % \
                (file_name, actual_file_md5sum, f.get('md5sum')))


def bam_check(metadata, seq_files):
    file_name_to_path = {}
    for f in seq_files:
        file_name_to_path[os.path.basename(f)] = f

    rg_info = {}  # RG ID as first key, value would be key-value of RG fields

    for f in metadata.get('files'):
        if not f.get('format') == 'BAM': continue

        file_name = f.get('name')
        bam_file = file_name_to_path.get(file_name)

        # run picard tool ValidateSamFile: b.12
        picard_validate_bam(bam_file)

        # get @RG lines from BAM header
        get_rg_from_bam_header(bam_file, rg_info)

    # perform checks on read group metadata
    rg_info_validation(metadata, rg_info)


def fastq_check(metadata, seq_files):
    return True


def run_validation(args):
    with open(args.metadata_json, 'r') as f:
        metadata = json.load(f)

    file_check(metadata, args.seq_files)

    bam_check(metadata, args.seq_files)

    fastq_check(metadata, args.seq_files)

    # reaching this point means no error found, reporting valid
    print(json.dumps({
        'valid': 'valid'
    }))


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-d", "--seq-files", dest="seq_files", help="Seq files to submit and process", type=str, nargs='+')
    parser.add_argument("-p", "--metadata-json", dest="metadata_json",
                        help="json file containing experiment, read_group and file information for sequence preprocessing")
    args = parser.parse_args()


    run_validation(args)
