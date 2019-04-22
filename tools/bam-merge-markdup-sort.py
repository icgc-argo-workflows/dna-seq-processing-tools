#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import subprocess
import argparse


def main():
    """ Main program """
    parser = argparse.ArgumentParser(description='Merge bam')
    parser.add_argument('tool', type=str, help='Tool to used',
                        choices=['samtools', 'picard-biobambam', 'sambamba', 'biobambam'])
    parser.add_argument('-i','--input-bams', dest='input_bams',
                        type=str, help='Input bam file', nargs='+', required=True)
    parser.add_argument('-o','--output-bam', dest='output_bam',
                        type=str, help='Output merged bam file', required=True)
    args = parser.parse_args()

    cmd = ''
    if args.tool == "samtools":
        merge = 'samtools merge - %s -f' % ' '.join(args.input_bams)
        sort = 'samtools sort - -o /dev/stdout'
        markdup = 'samtools markdup - %s' % (args.output_bam)
        cmd = '|'.join([merge, sort, markdup])

    elif args.tool == "picard-biobambam":
        merge = 'java -jar /tools/picard.jar MergeSamFiles I=%s O=/dev/stdout' % \
                    ' I='.join(args.input_bams)
        sort = 'java -jar /tools/picard.jar SortSam INPUT=/dev/stdin OUTPUT=/dev/stdout SORT_ORDER=coordinate'
        # Picard MarkDuplicates does not work with stdin as input, there seems to be a bug
        #markdup = 'java -jar /tools/picard.jar MarkDuplicates ASSUME_SORT_ORDER=coordinate I=/dev/stdin O=%s M=marked_dup_metrics.txt' % args.output_bam
        markdup = 'bammarkduplicates2 O=%s markthreads=16 M=duplicates_metrics.txt' % \
                    args.output_bam
        cmd = '|'.join([merge, sort, markdup])

    elif args.tool == "sambamba":
        merge = 'sambamba merge -t 16 -l 1 /dev/stdout %s' % ' '.join(args.input_bams)
        sort = 'sambamba sort -t 16 -l 1 -o /dev/stdout /dev/stdin'
        markdup = 'sambamba markdup -t 16 -l 5 /dev/stdin %s' % args.output_bam
        cmd = '|'.join([merge, sort, markdup])

    elif args.tool == "biobambam":
        raise NotImplementedError

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


if __name__ == "__main__":
    main()