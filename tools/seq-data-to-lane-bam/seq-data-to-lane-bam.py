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

import os
import subprocess
import sys
import json
import re
import hashlib
from argparse import ArgumentParser
from multiprocessing import cpu_count
import glob


def group_readgroup_by_filepair(seq_experiment_analysis):
    filepair_map_to_readgroup = {}

    # we assume information in seq_experiment_analysis has gone through
    # all validation checks in sequencing experiment submission
    # note: since advanced SONG validation is not ready, here we still validate uniqueness of
    #       read_group_id_in_bam and submitter_read_group_id
    read_group_id_in_bam_set = set()
    submitter_read_group_id_set = set()
    for rg in seq_experiment_analysis.get('read_groups'):
        rg['experiment'] = seq_experiment_analysis['experiment']  # let read group carry experiment
        rg['submitter_sample_id'] = seq_experiment_analysis['samples'][0]['submitterSampleId']  # let read group carry submitter_sample_id

        file_r1_r2 = (rg.get('file_r1'), rg.get('file_r2'))  # tuple

        if file_r1_r2 not in filepair_map_to_readgroup:
            filepair_map_to_readgroup[file_r1_r2] = {
                'format': 'BAM' if file_r1_r2[0].endswith('.bam') else 'FASTQ',
                'read_groups': [rg]
            }

        else:
            if filepair_map_to_readgroup[file_r1_r2]['format'] == 'FASTQ':  # shouldn't happen but just in case
                sys.exit('Error found: same FASTQ %s must not be associated to more than one read group\n' % \
                    ' and '.join(file_r1_r2) if file_r1_r2[1] else file_r1_r2[0])
            filepair_map_to_readgroup[file_r1_r2]['read_groups'].append(rg)

        # make sure no duplicate of read_group_id_in_bam (when populated)
        if rg.get('read_group_id_in_bam'):
            if rg['read_group_id_in_bam'] in read_group_id_in_bam_set:
                sys.exit('Error found: read_group_id_in_bam duplicated: %s' % rg['read_group_id_in_bam'])
            else:
                read_group_id_in_bam_set.add(rg['read_group_id_in_bam'])

        # make sure no duplicate of submitter_read_group_id
        if rg['submitter_read_group_id'] in submitter_read_group_id_set:
            sys.exit('Error found: submitter_read_group_id duplicated: %s' % rg['submitter_read_group_id'])
        else:
            submitter_read_group_id_set.add(rg['submitter_read_group_id'])

    return filepair_map_to_readgroup


def readgroup_id_to_fname(rg_id, input_bam_name='', study_id=None, donor_id=None, sample_id=None):
    friendly_rgid = "".join([ c if re.match(r"[a-zA-Z0-9\.\-_]", c) else "_" for c in rg_id ])
    # use original bam file name and rg_id to calculate the md5sum to avoid new lane bam file name collision
    # use white space (' ') to separate bam name and rg_id
    md5sum = hashlib.md5(("%s %s" % (input_bam_name, rg_id)).encode('utf-8')).hexdigest()

    if not sample_id or not donor_id or not study_id:
        sys.exit('Error: missing study/donor/sample ID in the provided metadata')

    return ".".join([study_id, donor_id, sample_id, friendly_rgid, md5sum, 'lane.bam'])


def bunzip2(fq_pair):
    bunzipped = []
    for f in fq_pair:
        if f and f.endswith('.bz2'):
            seq_file_bunzipped = os.path.join(
                                    os.environ.get("TMPDIR", ''),
                                    os.path.basename(re.sub(r'\.bz2$', '', f))
                                )
            cmd = 'bunzip2 -k -c %s > %s' % (f, seq_file_bunzipped)
            try:
                subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            except Exception as e:
                sys.exit("Unable to uncompress bz2 file: %s. Error: %s" % (f, e))

            bunzipped.append(seq_file_bunzipped)
        else:
            bunzipped.append(f)

    return bunzipped


def generate_ubam_from_fastq(fq_pair, readgroup, mem=None, study_id=None, donor_id=None, sample_id=None):
    jvm_Xmx = "-Xmx%sM" % mem if mem else ""
    fastq_pair = bunzip2(fq_pair)

    rg_args = [
        'READ_GROUP_NAME=%s' % readgroup['submitter_read_group_id'],
        'SAMPLE_NAME=%s' % readgroup['submitter_sample_id'],
        'LIBRARY_NAME=%s' % readgroup.get('library_name'),
        'PLATFORM_UNIT=%s' % readgroup['platform_unit'],
        'PLATFORM=%s' % readgroup['experiment'].get('platform')
    ]
    if readgroup.get('insert_size'):
        rg_args.append('PREDICTED_INSERT_SIZE=%s' % readgroup.get('insert_size'))
    if readgroup['experiment'].get('sequencing_center'):
        rg_args.append('SEQUENCING_CENTER=%s' % readgroup['experiment'].get('sequencing_center'))
    if readgroup['experiment'].get('platform_model'):
        rg_args.append('PLATFORM_MODEL=%s' % readgroup['experiment'].get('platform_model'))
    if readgroup['experiment'].get('sequencing_date'):
        rg_args.append('RUN_DATE=%s' % readgroup['experiment'].get('sequencing_date'))

    # convert paired end fastq to unaligned and lane level bam sorted by query name
    # convert readGroupId to filename friendly
    try:
        cmd = [
                'java', jvm_Xmx, '-jar', '/tools/picard.jar', 'FastqToSam', 'FASTQ=%s' % fastq_pair[0],
                'OUTPUT=%s' % readgroup_id_to_fname(readgroup['submitter_read_group_id'], '', study_id, donor_id, sample_id)
              ] + rg_args

        if fastq_pair[1]:
            cmd += ['FASTQ2=%s' % fastq_pair[1]]

        subprocess.run(cmd, check=True)

    except Exception as e:
        sys.exit("Error: %s. FastqToSam failed: %s\n" % (e,
                    ' and '.join(fastq_pair) if fastq_pair[1] else fastq_pair[0]))


def generate_ubams_from_bam(bam, readgroups, mem=None, study_id=None, donor_id=None, sample_id=None):
    # get input bam basename, remove extention to use as output subfolder name
    bam_base = os.path.splitext(os.path.basename(bam))[0]
    out_format = bam_base+'/%!.bam'
    os.mkdir(bam_base)
    cmd = ['samtools', 'split', '-f', '%s' % out_format, '-@ %s' % str(args.cpus), bam]
    try:
        subprocess.run(cmd, check=True)
    except Exception as e:
        sys.exit("Error: %s. 'samtools split' failed: %s\n" % (e, bam))

    # convert readGroupId to filename friendly
    # only process the lanes output for given input bam
    for bamfile in glob.glob(os.path.join(os.getcwd(), bam_base, "*.bam")):

        # remove file extension to get rg_id
        rg_id = os.path.splitext(os.path.basename(bamfile))[0]

        # let's make sure RG_ID in lane bam exists in readgroup metadata, either matching read_group_id_in_bam or submitter_read_group_id
        rg_id_found = False
        for rg in readgroups:
            if rg.get('read_group_id_in_bam') == rg_id or \
                    (not rg.get('read_group_id_in_bam') and rg['submitter_read_group_id'] == rg_id):
                rg_id_found = True
                break

        if not rg_id_found:
            sys.exit("Error: unable to find read group info for rg_id ('%s') in the supplied metadata (SONG Analysis)" % rg_id)

        os.rename(bamfile, os.path.join(os.getcwd(), readgroup_id_to_fname(rg_id, os.path.basename(bam), study_id, donor_id, sample_id)))


def filename_to_file(filenames: tuple, files: list) -> tuple:
    name_to_file_map = {}
    for f in files:
        name_to_file_map[os.path.basename(f)] = f

    return (name_to_file_map[filenames[0]], name_to_file_map[filenames[1]] if filenames[1] else None)


def main(args):
    with open(args.metadata_json, 'r') as f:
        seq_experiment_analysis = json.load(f)

    study_id = seq_experiment_analysis['studyId']
    donor_id = seq_experiment_analysis['samples'][0]['donor']['donorId']
    sample_id = seq_experiment_analysis['samples'][0]['sampleId']

    filepair_map_to_readgroup = group_readgroup_by_filepair(seq_experiment_analysis)

    for fp in filepair_map_to_readgroup:
        if filepair_map_to_readgroup[fp]['format'] == 'BAM':
            # for bam just need fp[0] since fp[1] is either the same as fp[0] or None
            generate_ubams_from_bam(filename_to_file(fp, args.seq_files)[0],
                                        filepair_map_to_readgroup[fp]['read_groups'],
                                        args.mem, study_id, donor_id, sample_id)
        else:
            generate_ubam_from_fastq(filename_to_file(fp, args.seq_files),
                                        filepair_map_to_readgroup[fp]['read_groups'][0],
                                        args.mem, study_id, donor_id, sample_id)  # FASTQ must be one read group


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-s", "--seq-files", dest="seq_files", required=True,
                        help="Seq files to process", type=str, nargs='+')
    parser.add_argument("-p", "--metadata-json", dest="metadata_json", required=True,
                        help="JSON file containing sequencing_experiment analysis")
    parser.add_argument("-d", "--max-discard-fraction",
                        dest="max_discard_fraction",
                        default=0.05, type=float,
                        help="Max fraction of reads allowed to be discarded when reverting aligned BAM to unaligned")
    parser.add_argument("-n", "--cpus", type=int, default=cpu_count())
    parser.add_argument("-m", "--mem", dest="mem", help="Maximal allocated memory in MB", type=int)
    args = parser.parse_args()

    main(args)
