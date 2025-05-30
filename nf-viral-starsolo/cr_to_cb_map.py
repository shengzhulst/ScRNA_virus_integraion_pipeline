#!/usr/bin/env python
# encoding: utf-8

import argparse
import pysam

arguments = argparse.ArgumentParser(description="Extract CR to CB map and UR to UB map")

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
    all_cr = set()
    all_cb = set()
    all_ur = set()
    all_ub = set()
    for read in samfile:
        if read.has_tag('CR'):
            cr = read.get_tag('CR')
            if read.has_tag('CB'):
                cb = read.get_tag('CB')
                if cr not in all_cr:
                    f.write(cr + '\t' + cb + '\n')
                    all_cr.add(cr)
                    all_cb.add(cb)
        if read.has_tag('UR'):
            ur = read.get_tag('UR')
            if read.has_tag('UB'):
                ub = read.get_tag('UB')
                if ur not in all_ur:
                    f1.write(ur + '\t' + ub + '\n')
                    all_ur.add(ur)
                    all_ub.add(ub)
print(f'Number of CR: {len(all_cr)}, number of CB: {len(all_cb)}')
print(f'Number of UR: {len(all_ur)}, number of UB: {len(all_ub)}')
samfile.close()

