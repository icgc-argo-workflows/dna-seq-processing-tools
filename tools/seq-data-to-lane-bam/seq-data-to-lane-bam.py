#!/usr/bin/env python3

import os
import subprocess
import sys
import json
import re
import glob
import datetime
from argparse import ArgumentParser

"""
Major steps:
- convert input Seq to unaligned BAM for each read group
"""

def main(args):
    with open(args.metadata_json, 'r') as f:
        metadata = json.load(f)
    picard = args.picard_jar
    cwd = os.getcwd()
    output = {}
    if args.seq_format == 'FASTQ':
        readGroups = metadata.get('read_groups')
        for rg in readGroups:
            readGroupId = rg.get('submitter_id')
            files = rg.get('files')
            file_with_path = []
            for _file in files:
                file_path = os.path.join(args.seq_files_dir, _file.get('name'))
                if not os.path.isfile(file_path): sys.exit('\n The file: %s do not exist!' % file_path)
                file_with_path.append(file_path)

            # detect whether there are more than two fastq files for each read_group
            if not len(file_with_path) == 2:
                sys.exit('\nThe number of fastq files is not equal to 2 for %s' % readGroupId)

            rg_args = ['READ_GROUP_NAME=%s' % readGroupId,
                       'SAMPLE_NAME=%s' % metadata.get('sample_submitter_id'),
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
                subprocess.run(['java', '-jar', picard,
                                'FastqToSam', 'FASTQ=%s' % file_with_path[0],
                                'FASTQ2=%s' % file_with_path[1],
                                'OUTPUT=%s' % os.path.join(cwd, rg_fname+'.lane.bam')] + rg_args, check=True)
            except Exception as e:
                sys.exit('\n%s: FastqToSam failed: %s and %s' % (e, file_with_path[0], file_with_path[1]))

    # the inputs are BAM
    elif args.seq_format == 'BAM':
        files = metadata.get('files')

        for _file in files:
            file_path = os.path.join(args.seq_files_dir, _file.get('name'))
            if not os.path.isfile(file_path): sys.exit('\n The file: %s do not exist!' % file_path)

            # Revert the bam to unaligned and lane level bam sorted by query name
            # Suggested options from: https://github.com/broadinstitute/picard/issues/849#issuecomment-313128088
            try:
                subprocess.run(['java', '-jar', picard,
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

    output['aligned_basename'] = '.'.join([metadata.get('sample_submitter_id'), str(metadata.get('read_group_count')), datetime.date.today().strftime("%Y%m%d"), 'wgs', 'grch38'])
    output['bundle_type'] = "lane_seq_submission"

    # write the parameter to stdout
    print(json.dumps(output))

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-i", "--seq_format", dest="seq_format",
                        help="Sequence format")
    parser.add_argument("-d", "--seq_files_dir", dest="seq_files_dir", help="Directory with seq files to submit and process")
    parser.add_argument("-p", "--metadata_json", dest="metadata_json",
                        help="json file containing experiment, read_group and file information for sequence preprocessing")
    parser.add_argument("-j", "--picard_jar", dest="picard_jar", default="/tools/picard.jar",
                        help="Picard jar file")
    args = parser.parse_args()

    main(args)