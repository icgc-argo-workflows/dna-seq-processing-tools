#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
  Copyright (C) 2021,  icgc-argo

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU Affero General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Affero General Public License for more details.

  You should have received a copy of the GNU Affero General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

  Authors:
    Linda Xiang
"""

import os
import sys
import subprocess
import json
import re
import hashlib
from argparse import ArgumentParser
from multiprocessing import cpu_count
import glob
import shutil
import csv

def run_cmd(cmd):
    proc = subprocess.Popen(
                cmd,
                shell=True,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE
            )
    stdout, stderr = proc.communicate()

    return (
        stdout.decode("utf-8").strip(),
        stderr.decode("utf-8").strip(),
        proc.returncode
    )

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

        # make sure no duplicate of read_group_id_in_bam (when populated) within the same bam
        if rg.get('read_group_id_in_bam'):
            bam_and_rg_id = '%s %s' % (rg.get('file_r1'), rg.get('read_group_id_in_bam'))
            if bam_and_rg_id in read_group_id_in_bam_set:
                sys.exit('Error found: read_group_id_in_bam duplicated in the same BAM: %s' % bam_and_rg_id)
            else:
                read_group_id_in_bam_set.add(bam_and_rg_id)

        # make sure no duplicate of submitter_read_group_id
        if rg['submitter_read_group_id'] in submitter_read_group_id_set:
            sys.exit('Error found: submitter_read_group_id duplicated: %s' % rg['submitter_read_group_id'])
        else:
            submitter_read_group_id_set.add(rg['submitter_read_group_id'])

    return filepair_map_to_readgroup

def filename_to_file(filenames: tuple, files: list) -> tuple:
    name_to_file_map = {}
    for f in files:
      name_to_file_map[os.path.basename(f)] = f

    return (name_to_file_map[filenames[0]], name_to_file_map[filenames[1]] if filenames[1] else None)

def readgroup_id_to_fname(rg_id, input_bam_name='', study_id=None, donor_id=None, sample_id=None):
    friendly_rgid = "".join([ c if re.match(r"[a-zA-Z0-9\.\-_]", c) else "_" for c in rg_id ])
    # use original bam file name and rg_id to calculate the md5sum to avoid new lane bam file name collision
    # use white space (' ') to separate bam name and rg_id
    md5sum = hashlib.md5(("%s %s" % (input_bam_name, rg_id)).encode('utf-8')).hexdigest()

    if not sample_id or not donor_id or not study_id:
        sys.exit('Error: missing study/donor/sample ID in the provided metadata')

    return ".".join([study_id, donor_id, sample_id, friendly_rgid, md5sum])


def generate_fastqs_from_bam(bam, readgroups, cpu=None, rgs_file_pair_map=dict(), study_id=None, donor_id=None, sample_id=None, temp_dir=None, out_dir=None):
    # get input bam basename, remove extention to use as output subfolder name
    bam_base = os.path.splitext(os.path.basename(bam))[0]
    out_format = bam_base+'/%!.bam'
    if os.path.exists(bam_base) and os.path.isdir(bam_base):
      shutil.rmtree(bam_base)
    os.mkdir(bam_base)
    cmd = ['samtools', 'split', '-f', '%s' % out_format, '-@ %s' % str(cpu), bam]

    stdout, stderr, returncode = run_cmd(" ".join(cmd))
    if returncode:
        sys.exit(f"Error: 'samtools split' failed.\nStdout: {stdout}\nStderr: {stderr}\n")

    # convert readGroupId to filename friendly
    # only process the lanes output for given input bam
    for lane_bam in glob.glob(os.path.join(os.getcwd(), bam_base, "*.bam")):

        # remove file extension to get rg_id
        rg_id = os.path.splitext(os.path.basename(lane_bam))[0]

        # let's make sure RG_ID in lane bam exists in readgroup metadata, either matching read_group_id_in_bam or submitter_read_group_id
        # otherwise, it should be a lane that's expected to be ignored
        rg_id_found = False
        for rg in readgroups:
            if rg.get('file_r1') == os.path.basename(bam) and (rg.get('read_group_id_in_bam') == rg_id or
                    (not rg.get('read_group_id_in_bam') and rg['submitter_read_group_id'] == rg_id)):
                rg_id_found = True
                # rgs_with_lane_bam_produced.add(rg['submitter_read_group_id'])
                break

        if rg_id_found:
            rg_id_fn = readgroup_id_to_fname(rg_id, os.path.basename(bam), study_id, donor_id, sample_id)
            if rg['is_paired_end']:
              cmd = f'samtools collate -uOn 128 {lane_bam} {temp_dir}/tmp_{rg_id_fn} | samtools fastq -O -F 0x900 -1 {out_dir}/{rg_id_fn}_R1.fq.gz -2 {out_dir}/{rg_id_fn}_R2.fq.gz -'
            else:
              cmd = f'samtools collate -uOn 128 {lane_bam} {temp_dir}/tmp_{rg_id_fn} | samtools fastq -O -F 0x900 -s {out_dir}/{rg_id_fn}_R1.fq.gz -'
            
            stdout, stderr, returncode = run_cmd(cmd)
            if returncode:
              sys.exit(f"Error: 'samtools collate and fastq' failed.\nStdout: {stdout}\nStderr: {stderr}\n")
            
            rgs_file_pair_map[rg['submitter_read_group_id']] = {
              'file_r1': os.path.join(os.getcwd(), f'{out_dir}/{rg_id_fn}_R1.fq.gz'),
              'file_r2': os.path.join(os.getcwd(), f'{out_dir}/{rg_id_fn}_R2.fq.gz') if os.path.isfile(f'{out_dir}/{rg_id_fn}_R2.fq.gz') else 'No_File'
            }

        else:  # ignore lane bam without read group information in metadata, just produce a warning message here
            print("WARNING: Ignore lane BAM '%s' (split from input BAM '%s') that has no corresponding read group in the metadata" %
                  (lane_bam, os.path.basename(bam)), file=sys.stderr)

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

def get_new_filename(fastq_old, rg_id_fn, r1_r2, outdir):
    if fastq_old.endswith('fq') or fastq_old.endswith('fastq'):
      fastq_new = os.path.join(os.getcwd(), outdir, f'{rg_id_fn}_{r1_r2}.fq')
    elif fastq_old.endswith('fq.gz') or fastq_old.endswith('fastq.gz'):
      fastq_new = os.path.join(os.getcwd(), outdir, f'{rg_id_fn}_{r1_r2}.fq.gz')
    else:
      sys.exit("Unsupported file format: %s." % fastq_old)
    
    return fastq_new
    


def main():
    """
    Python implementation of tool: seq-data-to-lane-fastq
    """

    parser = ArgumentParser(description='Tool: seq-data-to-lane-fastq')
    parser.add_argument("-s", "--seq-files", dest="seq_files", required=True,
                        help="Seq files to process", type=str, nargs='+')
    parser.add_argument("-p", "--metadata-json", dest="metadata_json", required=True,
                        help="JSON file containing sequencing_experiment analysis")
    parser.add_argument('-t', '--tempdir', dest='tempdir', type=str, default=".",
                        help='Specify directory for temporary files')
    parser.add_argument('-o', '--outdir', dest='outdir', type=str, default="out",
                        help='Specify directory for output files')
    parser.add_argument("-d", "--max-discard-fraction",
                        dest="max_discard_fraction",
                        default=0.05, type=float,
                        help="Max fraction of reads allowed to be discarded when reverting aligned BAM to unaligned")
    parser.add_argument("-n", "--cpus", type=int, default=cpu_count())

    args = parser.parse_args()

    if os.path.exists(args.outdir) and os.path.isdir(args.outdir):
      shutil.rmtree(args.outdir)
    os.mkdir(args.outdir)

    with open(args.metadata_json, 'r') as f:
      seq_experiment_analysis = json.load(f)

    study_id = seq_experiment_analysis['studyId']
    donor_id = seq_experiment_analysis['samples'][0]['donor']['donorId']
    sample_id = seq_experiment_analysis['samples'][0]['sampleId']

    filepair_map_to_readgroup = group_readgroup_by_filepair(seq_experiment_analysis)

    rgs_file_pair_map = dict()

    for fp in filepair_map_to_readgroup:
      if filepair_map_to_readgroup[fp]['format'] == 'BAM':
        # for bam just need fp[0] since fp[1] is either the same as fp[0] or None
        generate_fastqs_from_bam(filename_to_file(fp, args.seq_files)[0],
                                    filepair_map_to_readgroup[fp]['read_groups'],
                                    args.cpus, rgs_file_pair_map, study_id, donor_id, sample_id, args.tempdir, args.outdir)
      else: # FASTQ must be one read group
        fq_pair = filename_to_file(fp, args.seq_files)
        fastq_pair = bunzip2(fq_pair)
        rg = filepair_map_to_readgroup[fp]['read_groups'][0]
        rg_id = rg['submitter_read_group_id']
        rg_id_fn = readgroup_id_to_fname(rg_id, '', study_id, donor_id, sample_id)
        file_r1_new = get_new_filename(fastq_pair[0], rg_id_fn, "R1", args.outdir)
        os.symlink(os.path.abspath(fastq_pair[0]), file_r1_new)
        if fastq_pair[1]:
          file_r2_new = get_new_filename(fastq_pair[1], rg_id_fn, "R2", args.outdir)
          os.symlink(os.path.abspath(fastq_pair[1]), file_r2_new)
        else:
          file_r2_new = 'No_File'
        rgs_file_pair_map[rg_id] = {
          'file_r1': file_r1_new,
          'file_r2': file_r2_new
        }
        
    # now we check whether all read groups in metadata have produced lane bam
    rgs_missed_lane = set()
    for rg in seq_experiment_analysis['read_groups']:
        if rg['submitter_read_group_id'] not in rgs_file_pair_map:
            rgs_missed_lane.add(rg['submitter_read_group_id'])

    if rgs_missed_lane:  # throw error here if that happens
        sys.exit("Error: no lane BAM has been generated for some read groups: '%s'. "
                "Please make sure supplied sequencing files and metadata are correct." % "', '".join(rgs_missed_lane))

    with open(f'{args.outdir}/rgs_file_pair_map.csv', 'w', newline='') as f:
      csvwriter = csv.writer(f, delimiter=',')
      csvwriter.writerow(['read_group_id', 'file_r1', 'file_r2'])
      for k,v in rgs_file_pair_map.items():
        csvwriter.writerow([k, v['file_r1'], v['file_r2']])



if __name__ == "__main__":
    main()

