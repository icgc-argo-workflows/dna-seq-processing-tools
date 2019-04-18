#!/usr/bin/env python
# -*- coding: utf-8 -*-

import subprocess
import argparse

def main():
    """ Main program """
    parser = argparse.ArgumentParser(description='BWA alignment')
    parser.add_argument('-i','--input-bam', dest='input_bam', type=str, help='Input bam file', required=True)
    parser.add_argument('-r','--ref-genome', dest='ref_genome', type=str, help='Fa.gz reference genome file', required=True)
    parser.add_argument('-o','--output', dest='output', type=str, help='Output bam file', required=True)
    args = parser.parse_args()

    bam_header = subprocess.Popen(['samtools','view','-H',args.input_bam],stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()[0].strip()
    rg_header = []
    for line in bam_header.split('\n'):
        if "@RG" in line: rg_header.append(line)


    bam_to_fastq = 'samtools fastq %s' % args.input_bam

    #Command with header
    #alignment = 'bwa mem -R %s %s -' % (''.join(rg_header), args.ref_genome)

    #Command without header
    alignment = 'bwa mem %s -' % (args.ref_genome)

    subprocess.call(' '.join([bam_to_fastq,'|',alignment,'>',args.output]),shell=True)

if __name__ == "__main__":
    main()