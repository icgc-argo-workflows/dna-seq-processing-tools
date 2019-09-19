#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import subprocess
import argparse
from multiprocessing import cpu_count
import json


def main():
    """ Main program """
    parser = argparse.ArgumentParser(description='Merge and markdup')
    parser.add_argument('-i','--input-bams', dest='input_bams',
                        type=str, help='Input bam file', nargs='+', required=True)
    parser.add_argument('-b','--output-base', dest='output_base',
                        type=str, help='Output merged file basename', required=True)
    parser.add_argument('-r', '--reference', dest='reference',
                        type=str, help='reference fasta', required=True)
    parser.add_argument("-n", "--cpus", dest='cpus', type=int, default=cpu_count())
    parser.add_argument("-d", "--mdup", dest='mdup', action='store_true')
    parser.add_argument("-o", "--output-format", dest='output_format', nargs='+', default=['cram'])
    args = parser.parse_args()

    cmd = []

    if args.mdup:
        if "bam" in args.output_format and "cram" in args.output_format:
            markdup = 'bammarkduplicates2 markthreads=%s O=/dev/stdout M=%s indexfilename=%s index=1 I=%s ' % \
                      (str(args.cpus), args.output_base + ".duplicates-metrics.txt", args.output_base + ".bam.bai",
                       ' I='.join(args.input_bams))
            tee = 'tee %s ' % (args.output_base+".bam")
            cram = 'samtools view -C -T %s -@ %s /dev/stdin | tee %s ' % (
            args.reference, args.cpus, args.output_base + ".cram")
            crai = 'samtools index -@ %s /dev/stdin %s' % (args.cpus, args.output_base + ".cram.crai")
            cmd.append('|'.join([markdup, tee, cram, crai]))
        elif "bam" in args.output_format and not "cram" in args.output_format:
            markdup = 'bammarkduplicates2 markthreads=%s O=%s M=%s indexfilename=%s  index=1 I=%s ' % \
                      (str(args.cpus), args.output_base + ".bam", args.output_base + ".duplicates-metrics.txt",
                       args.output_base + ".bam.bai", ' I='.join(args.input_bams))
            cmd.append(markdup)
        elif not "bam" in args.output_format and "cram" in args.output_format:
            markdup = 'bammarkduplicates2 markthreads=%s O=/dev/stdout M=%s I=%s ' % \
                      (str(args.cpus), args.output_base + ".duplicates-metrics.txt", ' I='.join(args.input_bams))
            cram = 'samtools view -C -T %s -@ %s /dev/stdin ' % (args.reference, args.cpus)
            tee = 'tee %s ' % (args.output_base + ".cram")
            crai = 'samtools index -@ %s /dev/stdin %s' % (args.cpus, args.output_base + ".cram.crai")
            cmd.append('|'.join([markdup, cram, tee, crai]))

    elif not args.mdup:
        if "bam" in args.output_format and "cram" in args.output_format:
            merge = 'samtools merge -uf -@ %s /dev/stdout %s | tee %s' % (args.cpus, ' '.join(args.input_bams), args.output_base + ".bam")
            cram = 'samtools view -C -T %s -@ %s /dev/stdin | tee %s ' % (args.reference, args.cpus, args.output_base + ".cram")
            crai = 'samtools index -@ %s /dev/stdin %s' % (args.cpus, args.output_base + ".cram.crai")
            cmd.append('|'.join([merge, cram, crai]))
            bai = 'samtools index -@ %s %s %s ' % (args.cpus, args.output_base + ".bam", args.output_base + ".bam.bai")
            cmd.append(bai)
        elif "bam" in args.output_format and not "cram" in args.output_format:
            merge = 'samtools merge -uf -@ %s /dev/stdout %s | tee %s' % (args.cpus, ' '.join(args.input_bams), args.output_base + ".bam")
            bai = 'samtools index -@ %s /dev/stdin %s' % (args.cpus, args.output_base + ".bam.bai")
            cmd.append('|'.join([merge, bai]))
        elif not "bam" in args.output_format and "cram" in args.output_format:
            merge = 'samtools merge -uf -@ %s /dev/stdout %s ' % (args.cpus, ' '.join(args.input_bams))
            cram = 'samtools view -C -T %s -@ %s /dev/stdin | tee %s ' % (args.reference, args.cpus, args.output_base + ".cram")
            crai = 'samtools index -@ %s /dev/stdin %s' % (args.cpus, args.output_base + ".cram.crai")
            cmd.append('|'.join([merge, cram, crai]))


    for c in cmd:
        stdout, stderr, p, success = '', '', None, True
        try:
            p = subprocess.Popen([c],
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 shell=True)
            p.communicate()
        except Exception as e:
            print('Execution failed: %s' % e, file=sys.stderr)
            success = False

        if p and p.returncode != 0:
            print('Execution failed, none zero code returned.', file=sys.stderr)
            success = False


        if not success:
            sys.exit(p.returncode if p.returncode else 1)

    # write the parameter to stdout
    output = {"bundle_type": "dna_alignment"}
    print(json.dumps(output))


if __name__ == "__main__":
    main()
