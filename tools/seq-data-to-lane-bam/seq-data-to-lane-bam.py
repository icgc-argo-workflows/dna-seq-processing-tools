#!/usr/bin/env python3

import os
import subprocess
import sys
import json
import re
import uuid
import glob
import datetime
from argparse import ArgumentParser


def get_uuid5(bid, fid):
    uuid5 = str(uuid.uuid5(uuid.UUID("6ba7b810-9dad-11d1-80b4-00c04fd430c8"), "%s/%s" % (bid, fid)))
    return uuid5

def run_cmd(cmd):
    p, success = None, True
    try:
        p = subprocess.run(cmd,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             shell=True)
    except Exception as e:
        print('Execution failed: %s' % e)
        success = False

    if p and p.returncode != 0:
        print('Error occurred, return code: %s. Details: %s' % \
                (p.returncode, p.stderr.decode("utf-8")), file=sys.stderr)
        success = False

    if not success:
        sys.exit(p.returncode if p.returncode else 1)

    return


def main(args):
    with open(args.metadata_json, 'r') as f:
        metadata = json.load(f)
    cwd = os.getcwd()
    output = {}
    if metadata.get("input_seq_format") == 'FASTQ':
        readGroups = metadata.get('read_groups')
        for rg in readGroups:
            readGroupId = rg.get('submitter_read_group_id')
            files = rg.get('files')
            file_with_path = []
            for _file in files:
                for seq_file in args.seq_files:
                    if _file.get('name') != os.path.basename(seq_file): continue
                    if not os.path.isfile(seq_file): sys.exit('\n The file: %s do not exist!' % seq_file)
                    if seq_file.endswith(".bz2"):
                        seq_file_unzip = os.path.join(os.environ["TMPDIR"], _file.get('name').replace('.bz2', ''))
                        cmd = 'bunzip2 -k -c %s > %s' % (seq_file, seq_file_unzip)
                        run_cmd(cmd)
                        file_with_path.append(seq_file_unzip)
                    else:
                        file_with_path.append(seq_file)

            # detect whether there are more than two fastq files for each read_group
            if not len(file_with_path) == 2:
                sys.exit('\nThe number of fastq files is not equal to 2 for %s' % readGroupId)

            rg_args = ['READ_GROUP_NAME=%s' % readGroupId,
                       'SAMPLE_NAME=%s' % metadata.get('submitter_sample_id'),
                       'LIBRARY_NAME=%s' % rg.get('library_name'),
                       'PLATFORM_UNIT=%s' % rg.get('platform_unit'),
                       'PLATFORM=%s' % metadata.get('platform')]
            if metadata.get('sequencing_center') and str(metadata.get('sequencing_center')) != '':
                rg_args.append('SEQUENCING_CENTER=%s' % metadata.get('sequencing_center'))
            if rg.get('insert_size') and isinstance(rg.get('insert_size'), int):
                rg_args.append('PREDICTED_INSERT_SIZE=%s' % rg.get('insert_size'))
            if metadata.get('platform_model') and str(metadata.get('platform_model')) != '':
                rg_args.append('PLATFORM_MODEL=%s' % metadata.get('platform_model'))


            # convert pair end fastq to unaligned and lane level bam sorted by query name
            # convert readGroupId to filename friendly
            rg_fname = "".join([ c if re.match(r"[a-zA-Z0-9\-_]", c) else "_" for c in readGroupId ])
            try:
                subprocess.run(['java', '-jar', '/tools/picard.jar',
                                'FastqToSam', 'FASTQ=%s' % file_with_path[0],
                                'FASTQ2=%s' % file_with_path[1],
                                'OUTPUT=%s' % os.path.join(cwd, rg_fname+'.lane.bam')] + rg_args, check=True)
            except Exception as e:
                sys.exit('\n%s: FastqToSam failed: %s and %s' % (e, file_with_path[0], file_with_path[1]))

    # the inputs are BAM
    elif metadata.get("input_seq_format") == 'BAM':
        files = metadata.get('files')

        for _file in files:
            for seq_file in args.seq_files:
                if _file.get('name') != os.path.basename(seq_file): continue
                file_path = seq_file
            if not os.path.isfile(file_path): sys.exit('\n The file: %s do not exist!' % file_path)

            # Revert the bam to unaligned and lane level bam sorted by query name
            # Suggested options from: https://github.com/broadinstitute/picard/issues/849#issuecomment-313128088
            try:
                subprocess.run(['java', '-jar', '/tools/picard.jar',
                                'RevertSam',
                                'I=%s' % file_path,
                                'SANITIZE=true',
                                'ATTRIBUTE_TO_CLEAR=XT',
                                'ATTRIBUTE_TO_CLEAR=XN',
                                'ATTRIBUTE_TO_CLEAR=AS',
                                'ATTRIBUTE_TO_CLEAR=OC',
                                'ATTRIBUTE_TO_CLEAR=OP',
                                'SORT_ORDER=queryname',
                                'RESTORE_ORIGINAL_QUALITIES=true',
                                'REMOVE_DUPLICATE_INFORMATION=true',
                                'REMOVE_ALIGNMENT_INFORMATION=true',
                                'OUTPUT_BY_READGROUP=true',
                                'VALIDATION_STRINGENCY=LENIENT',
                                'MAX_DISCARD_FRACTION=%s' % args.max_discard_fraction,
                                'O=%s' % cwd], check=True)
            except Exception as e:
                sys.exit('\n%s: RevertSam failed: %s' %(e, file_path))

            for filename in glob.glob(os.path.join(cwd, "*.bam")):
                # convert readGroupId to filename friendly
                readGroupId = os.path.basename(filename).replace(".bam", "")
                rg_fname = "".join([c if re.match(r"[a-zA-Z0-9\-_]", c) else "_" for c in readGroupId])
                os.rename(filename, os.path.join(cwd, rg_fname+".lane.bam"))

    else:
        sys.exit('\n%s: Input files format are not FASTQ or BAM')

    uuid_prefix = get_uuid5(metadata.get('program_id'), metadata.get('submitter_sample_id'))
    output['aligned_basename'] = '.'.join([uuid_prefix, str(metadata.get('read_group_count')), datetime.date.today().strftime("%Y%m%d"), 'wgs', 'grch38'])
    output['bundle_type'] = "lane_seq_submission"

    # write the parameter to stdout
    print(json.dumps(output))

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-d", "--seq_files", dest="seq_files", help="Seq files to submit and process", type=str, nargs='+')
    parser.add_argument("-p", "--metadata_json", dest="metadata_json",
                        help="json file containing experiment, read_group and file information for sequence preprocessing")
    parser.add_argument("-m", "--max_discard_fraction",
                        dest="max_discard_fraction",
                        default=0.05, type=float,
                        help="Max fraction of reads allowed to be discarded when reverting aligned BAM to unaligned")
    args = parser.parse_args()

    main(args)
