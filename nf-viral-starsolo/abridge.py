#!/usr/bin/env python
# encoding: utf-8

import pandas as pd
import argparse

min_reads = 2

arguments = argparse.ArgumentParser(description="Remove insertions with fewer than min read supports")

arguments.add_argument(
    "-i",
    required=True,
    type=str,
    help="full tsv file",
)

arguments.add_argument(
    "-o",
    required=True,
    type=str,
    help="output tsv file",
)

arguments.add_argument(
    "-a",
    required=True,
    type=str,
    help="abridged tsv file",
)

args = arguments.parse_args()

output_file = args.o
abridged_file = args.a
input_tsv = args.i

# write abridged tsv
df = pd.read_csv(input_tsv, sep="\t")

#df = df[ (df.hits <= max_hits) & (df.total >= min_reads)]
df = df[ df.total >= min_reads ]
df.to_csv(output_file, sep="\t", index=False)

readnames = df.pop('readnames')
df = pd.concat([df, readnames], axis=1)
df.to_csv(abridged_file, sep="\t", index=False)
