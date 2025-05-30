#!/usr/bin/env python
# encoding: utf-8

import argparse
import pysam

arguments = argparse.ArgumentParser(description="Extract read name to barcode and UMI map")

arguments.add_argument(
    "-b",
    required=True,
    type=str,
    help="BAM file",
)
arguments.add_argument(
    "-o",
    required=True,
    type=str,
    help="output barcode file",
)
arguments.add_argument(
    "-u",
    required=True,
    type=str,
    help="output UMI file",
)

args = arguments.parse_args()

bam = args.b
output_file = args.o
umi_file = args.u

samfile = pysam.AlignmentFile(bam, 'rb')
with open(output_file, 'w') as f, open(umi_file, 'w') as f1:
    all_read1 = set()
    all_read2 = set()
    for read in samfile:
        qname = read.query_name
        if read.has_tag('CR'):
            cr = read.get_tag('CR')
            if qname not in all_read1:
                f.write(qname + '\t' + cr + '\n')
                all_read1.add(qname)
        if read.has_tag('UR'):
            ur = read.get_tag('UR')
            if qname not in all_read2:
                f1.write(qname + '\t' + ur + '\n')
                all_read2.add(qname)
samfile.close()

