#!/usr/bin/env python
import re
import sys
from glob import glob
import numpy as np
import pandas as pd
from sam_dsbra_grapher_ling import *

to_scan_dir = sys.argv[1]
summary_fns = []
for root, dirs, files in os.walk(to_scan_dir):
    for f in files:
        if "summary_table.txt" in f:
            summary_fns.append(os.path.join(root, f))
samples = parse_files(summary_fns)

df_out = pd.DataFrame(columns=["experiment_name", "mutation_count_cutoff", 
                              "analyzed_seqs", "wildtype", "mismatch", 
                              "insertion", "deletion", "compound", 
                              "median_deletion_size", "median_insertion_size",
                              "median_repair_size"])
#analyzed_seqs = []
#mutation_count_cutoffs = []
#num_deletion = []
#num_insertion = []
#num_compound = []
#num_wt = []
#median_deletion_size = []
#median_insertion_size = []
#median_repair_size = []

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

    # median insertion size
    dft = insertion_len_dist(df)
    median_ins_size = dft["insertion_length"].median()

    # median repair size
    dft = repair_len_dist(df)
    median_repair_size = dft["repair_size"].median()

    # experiment name
    exp_name = re.findall("(.*)_first100k", summary_name)[0]

    row = [exp_name, count_cutoff, len(df)]\
          + list(dft_mut_freq["count"].values)\
          + [median_del_size, median_ins_size, median_repair_size]

    df_out.loc[idx] = row