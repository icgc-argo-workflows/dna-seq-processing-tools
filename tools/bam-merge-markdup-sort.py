#!/usr/bin/env python
# -*- coding: utf-8 -*-

import subprocess
import argparse
import psutil

def main():
    """ Main program """
    parser = argparse.ArgumentParser(description='Merge bam')
    parser.add_argument('tool', type=str, help='Tool to used', choices=['samtools','picard','sambamba','biobambam'])
    parser.add_argument('-i','--input-bam', dest='input_bam', type=str, help='Input bam file', nargs='+',required=True)
    parser.add_argument('-o','--output-bam', dest='output_bam', type=str, help='Output merged bam file', required=True)
    args = parser.parse_args()

    if args.tool == "samtools":
        merge = 'samtools merge /dev/stdout %s -f' % (' '.join(args.input_bam))
        sort = 'samtools sort -o /dev/stdout -'
        markdup = 'samtools markdup - %s' % (args.output_bam)
        subprocess.call('|'.join([merge,sort,markdup]), shell=True)
    if args.tool == "picard":
        raise NotImplementedError
    if args.tool == "sambamba":
        raise NotImplementedError
    if args.tool == "biobambam":
        raise NotImplementedError
    return


if __name__ == "__main__":
    main()