#!/usr/bin/env python
# encoding: utf-8

import argparse
import os
import numpy as np

import pandas as pd

arguments = argparse.ArgumentParser(description="Add cellranger corrected cell barcodes")

arguments.add_argument(
    "-i",
    required=True,
    type=str,
    help="summary tsv file",
)
arguments.add_argument(
    "-o",
    required=True,
    type=str,
    help="output tsv file",
)
arguments.add_argument(
    "-u",
    required=True,
    type=str,
    help="output with unique UMI",
)
arguments.add_argument(
    "-f",
    required=True,
    type=str,
    help="output filtered to barcode",
)
arguments.add_argument(
    "--cr_map",
    required=True,
    type=str,
    help="CR to CB map file",
)
arguments.add_argument(
    "--cr_map_vif",
    required=True,
    type=str,
    help="CR to CB map file from starsolo",
)
arguments.add_argument(
    "--id_cr_map",
    required=True,
    type=str,
    help="id to CR map file",
)
arguments.add_argument(
    "--ur_map",
    required=True,
    type=str,
    help="UR to UB map file",
)
arguments.add_argument(
    "--ur_map_vif",
    required=True,
    type=str,
    help="UR to UB map file from starsolo",
)
arguments.add_argument(
    "--id_ur_map",
    required=True,
    type=str,
    help="id to UR map file",
)
arguments.add_argument(
    "--barcode",
    required=True,
    type=str,
    help="filtered barcode list",
)

args = arguments.parse_args()

output_file = args.o
output_uniq = args.u
output_filtered = args.f
input_tsv = args.i
cr_map = args.cr_map
cr_map_vif = args.cr_map_vif
id_cr_map = args.id_cr_map
ur_map = args.ur_map
ur_map_vif = args.ur_map_vif
id_ur_map = args.id_ur_map
fb = args.barcode


input_tsv = pd.read_csv(input_tsv, sep="\t")
input_tsv = input_tsv.assign(readnames=input_tsv['readnames'].str.split(',')).explode('readnames')
readnames = input_tsv['readnames'].to_numpy()

fb = np.array(pd.read_csv(fb, header=None)[0])

cr_to_cb_map = {}
with open(cr_map) as f:
    for line in f:
        tmp = line.rstrip().split('\t')
        cr_to_cb_map[tmp[0]] = tmp[1]

cr_to_cb_map_vif = {}
with open(cr_map_vif) as f:
    for line in f:
        tmp = line.rstrip().split('\t')
        cr_to_cb_map_vif[tmp[0]] = tmp[1]

id_to_cr_map = {}
with open(id_cr_map) as f:
    for line in f:
        tmp = line.rstrip().split('\t')
        id_to_cr_map[tmp[0]] = tmp[1]

cr_list = []
count_no_map = 0
for x in readnames:
    if x in id_to_cr_map:
        cr_list.append(id_to_cr_map[x])
    else:
        cr_list.append('')
        count_no_map += 1
print(f'{count_no_map} ids are not in the id to CR map')

cb_list = []
count_no_map = 0
count_no_map_vif = 0
for x in cr_list:
    if x in cr_to_cb_map:
        cb_list.append(cr_to_cb_map[x])
    else:
        if x == '':
            cb_list.append(x + '-1')
        else:
            if cr_to_cb_map_vif[x] == '-':
                cb_list.append(x + '-1')
                count_no_map_vif += 1
            else:
                cb_list.append(cr_to_cb_map_vif[x])
            count_no_map += 1
print(f'{count_no_map} CRs are not in the CR to CB map, out of which, {count_no_map_vif} CRs are not in the starsolo CR to CB map either')

ur_to_ub_map = {}
with open(ur_map) as f:
    for line in f:
        tmp = line.rstrip().split('\t')
        ur_to_ub_map[tmp[0]] = tmp[1]

ur_to_ub_map_vif = {}
with open(ur_map_vif) as f:
    for line in f:
        tmp = line.rstrip().split('\t')
        ur_to_ub_map_vif[tmp[0]] = tmp[1]

id_to_ur_map = {}
with open(id_ur_map) as f:
    for line in f:
        tmp = line.rstrip().split('\t')
        id_to_ur_map[tmp[0]] = tmp[1]

ur_list = []
count_no_map = 0
for x in readnames:
    if x in id_to_ur_map:
        ur_list.append(id_to_ur_map[x])
    else:
        ur_list.append('')
        count_no_map += 1
print(f'{count_no_map} ids are not in the id to UR map')

ub_list = []
count_no_map = 0
count_no_map_vif = 0
for x in ur_list:
    if x in ur_to_ub_map:
        ub_list.append(ur_to_ub_map[x])
    else:
        if x == '':
            ub_list.append(x + '-1')
        else:
            if ur_to_ub_map_vif[x] == '-':
                ub_list.append(x + '-1')
                count_no_map_vif += 1
            else:
                ub_list.append(ur_to_ub_map_vif[x])
            count_no_map += 1
print(f'{count_no_map} URs are not in the UR to UB map, out of which, {count_no_map_vif} URs are not in the starsolo UR to UB map either')

input_tsv.insert(0, 'cellranger_barcode', cb_list)
input_tsv.insert(1, 'cellranger_UMI', ub_list)
input_tsv = input_tsv.sort_values(by=['cellranger_barcode', 'cellranger_UMI'])
input_tsv = input_tsv[(input_tsv['cellranger_barcode'] != '') & (input_tsv['cellranger_UMI'] != '')]
input_tsv.to_csv(output_file, sep="\t", index=False)

input_tsv = input_tsv[~input_tsv.duplicated(subset=['cellranger_barcode', 'cellranger_UMI', 'contig'], keep=False)]
input_tsv.to_csv(output_uniq, sep='\t', index=False)

input_tsv = input_tsv[np.isin(input_tsv['cellranger_barcode'].to_numpy(), fb)]
input_tsv.to_csv(output_filtered, sep="\t", index=False)

