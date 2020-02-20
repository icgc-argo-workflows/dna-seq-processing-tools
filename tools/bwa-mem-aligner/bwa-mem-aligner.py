#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess
import argparse
from multiprocessing import cpu_count
import sys
import os
import json
import re


def get_read_group_info(metadata_file, input_bam):
    if metadata_file and metadata_file != 'NO_FILE':
        with open(metadata_file, 'r') as f:
            metadata = json.load(f)

        # input bam naming convention: <read_group_id>.<md5sum>.lane.bam
        matched = re.match(r'^(.+?)\.[a-f0-9]{32}\.lane\.bam$', input_bam)
        if matched:
            read_group_id = matched.group(1)
        else:
            sys.exit("Error: input bam does not follow naming convention '<read_group_id>.<md5sum>.lane.bam': %s" % input_bam)
    else:
        return {}

    read_group = {}
    for rg in metadata['read_groups']:
        if rg['submitter_read_group_id'] == read_group_id:
            read_group = rg
            break
    if not read_group:
        sys.exit("Error: unable to find read group info for '%s' in the supplied metadata" % read_group_id)

    experiment = metadata['experiment']

    read_group_info = {
        'ID': read_group_id,
        'SM': metadata['samples'][0]['sampleId'],
        'LB': read_group['library_name'],
        'PU': read_group['platform_unit']
    }

    if read_group.get('insert_size'):
        read_group_info.update({'PI': read_group['insert_size']})
    if read_group.get('sample_barcode'):
        read_group_info.update({'BC': read_group['sample_barcode']})
    if experiment.get('sequencing_center'):
        read_group_info.update({'CN': experiment['sequencing_center']})
    if experiment.get('platform'):
        read_group_info.update({'PL': experiment['platform']})
    if experiment.get('platform_model'):
        read_group_info.update({'PM': experiment['platform_model']})
    if experiment.get('sequencing_date'):
        read_group_info.update({'DT': experiment['sequencing_date']})

    description = '|'.join([
                                experiment['library_strategy'],
                                metadata['studyId'],
                                metadata['samples'][0]['specimenId'],
                                metadata['samples'][0]['donor']['donorId'],
                                metadata['samples'][0]['specimen']['specimenType'],
                                metadata['samples'][0]['specimen']['tumourNormalDesignation']
                            ])

    read_group_info.update({'DS': description})

    return read_group_info


def main():
    """ Main program """
    parser = argparse.ArgumentParser(description='BWA alignment')
    parser.add_argument('-i','--input-bam', dest='input_bam', type=str,
                        help='Input bam file', required=True)
    parser.add_argument('-r','--ref-genome', dest='ref_genome', type=str,
                        help='Reference genome file (eg, .fa.gz), make sure BWA index files '
                             '(eg. .alt, .ann, .bwt etc) are all present at the same location', required=True)
    parser.add_argument('-o','--aligned_lane_prefix', dest='aligned_lane_prefix', type=str,
                        help='Output aligned lane bam file prefix', required=True)
    parser.add_argument("-n", "--cpus", type=int, default=cpu_count())
    parser.add_argument('-m','--metadata', dest='metadata', type=str,
                        help='Sequencing experiment metadata')
    args = parser.parse_args()

    # retrieve read_group_info from metadata
    read_group_info = get_read_group_info(args.metadata, args.input_bam)

    # retrieve the @RG from BAM header
    try:
        header = subprocess.check_output(['samtools', 'view', '-H', args.input_bam])

    except Exception as e:
        sys.exit('\n%s: Retrieve BAM header failed: %s' % (e, args.input_bam))

    # get @RG
    header_array = header.decode('utf-8').rstrip().split('\n')
    rg_array = []
    for line in header_array:
        if not line.startswith("@RG"): continue
        rg_array.append(line.rstrip().replace('\t', '\\t'))

    if not len(rg_array) == 1: sys.exit('\n%s: The input bam should only contain one readgroup ID: %s' % args.input_bam)

    if read_group_info:  # use what's in metadata instead of in BAM header
        rg_kv = [ '@RG' ] + [ '%s:%s' % (k, v) for k, v in read_group_info.items() ]
        rg_array = [ '\\t'.join(rg_kv) ]

    sort_qname = 'samtools sort -l 0 -n -O BAM -@ %s %s ' % (str(args.cpus), args.input_bam)

    # only output paired reads, discarding singletons, supplementary and secondary reads.
    bam2fastq = ' samtools fastq -O -F 0x900 -0 /dev/null -s /dev/null -@ %s - ' % (str(args.cpus))

    #Command with header
    alignment = ' bwa mem -K 100000000 -Y -t %s -p -R "%s" %s - ' % (str(args.cpus), rg_array[0], args.ref_genome)

    # Sort the SAM output by coordinate from bwa and save to BAM file
    sort_coordinate = ' samtools sort --no-PG -O BAM -@ %s -o %s /dev/stdin' % (str(args.cpus), args.aligned_lane_prefix + "." + os.path.basename(args.input_bam))

    cmd = ' | '.join([sort_qname, bam2fastq, alignment, sort_coordinate])

    try:
        subprocess.run([cmd], shell=True, check=True)

    except Exception as e:
        sys.exit('\nExecution failed: %s' % e)


if __name__ == "__main__":
    main()
