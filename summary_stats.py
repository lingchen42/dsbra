#!/usr/bin/env python
import re
import sys
from glob import glob
import numpy as np
import pandas as pd
from sam_dsbra_grapher_ling import *

gene_list = ["NFR", "CAN1Y", "UBP4"]

to_scan_dir = sys.argv[1]
summary_fns = []
for root, dirs, files in os.walk(to_scan_dir):
    for f in files:
        if "summary_table.txt" in f:
            summary_fns.append(os.path.join(root, f))
samples = parse_files(summary_fns)

df_out = pd.DataFrame(columns=["experiment_name", "gene", "count_cutoff", 
                              "analyzed_seqs", "wildtype", "mismatch", 
                              "insertion", "deletion", "compound", 
                              "median_deletion_size", "mean_deletion_size",
                              "median_insertion_size", "mean_insertion_size",
                              "median_repair_size", "mean_repair_size"])

for idx, summary_fn in enumerate(summary_fns):
    summary_name = os.path.basename(summary_fn).replace('.txt', '')
    cols = [col for col in samples[summary_fn][0].keys() if col]
    df = pd.DataFrame.from_records(samples[summary_fn], columns=cols) 

    # count cutoff
    print("Using poisson distriution of 300 colonies, p=0.001 as cutoff...")
    p, num_colony = 0.001, 300
    df, count_cutoff = apply_count_cutoff(df, p, num_colony)

    # mutation type frequency
    dft_mut_freq = get_mutation_event_freq(df, summary_fn)

    # median deletion size
    dft = del_len_dist(df)
    median_del_size = dft["deletion_length"].median()
    mean_del_size = dft["deletion_length"].mean()

    # median insertion size
    dft = insertion_len_dist(df)
    median_ins_size = dft["insertion_length"].median()
    mean_ins_size = dft["insertion_length"].mean()

    # median repair size
    dft = repair_len_dist(df)
    median_repair_size = dft["repair_size"].median()
    mean_repair_size = dft["repair_size"].mean()

    # experiment name
    exp_name = re.findall("(.*)_first100k", summary_name)[0]

    # gene
    eles = summary_fn.split("/")
    for ele in eles:
        if ele.upper() in gene_list:
            gene = ele

    row = [exp_name, gene, count_cutoff, df["count"].sum()]\
          + list(dft_mut_freq["count"].values)\
          + [median_del_size, mean_del_size, 
            median_ins_size, mean_ins_size,
            median_repair_size, median_repair_size]

    df_out.loc[idx] = row

df_out["insertion_frac"] = df_out["insertion"] / df_out["analyzed_seqs"]
df_out["deletion_frac"] = df_out["deletion"] / df_out["analyzed_seqs"]

outfn = os.path.join(to_scan_dir, "summary_stats.csv")
df_out.to_csv(outfn)
print("Complete. Results are written in ", outfn)