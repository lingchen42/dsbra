#!/usr/bin/env python
'''
This script is to visualize the mutation events in the output of sam_dsbra_interface.py
Usage:
    sam_dsbra_grapher_ling.py -h

Ling Chen 2018-05-29
'''
from __future__ import division
import csv
import ast
import argparse
from datetime import datetime
import pandas as pd
import os
import subprocess
import itertools
import re


def parse_files(summary_fns):
    '''
    Parse mutation summary files.
    Args:
        summary_fns: the output txt files from sam_dsbra_interface.py
    Return:
        samples, sample_count_list
    '''
    #samples: This dictionary stores the data in the output txt file
    samples = {}
    #sample_count_list: This dictionary stores the number analyzed sequences
    #                   in each file.
    sample_count_list = {}

    for summary_fn in summary_fns:
        with open(summary_fn,'r') as sf:
            sample_count_list[summary_fn] = int(sf.readline()[1:])

            data = list(csv.DictReader(sf,delimiter='\t'))
            #formatting of all fields appropriately
            for x in data:
                x['insertion_seqs'] = ast.literal_eval(x['insertion_seqs'])
                x['deletion_lens'] = ast.literal_eval(x['deletion_lens'])
                x['num_transitions'] = int(x['num_transitions'])
                x['num_deletions'] = int(x['num_deletions'])
                x['num_insertions'] = int(x['num_insertions'])
                x['new_mismatch_bases'] = ast.literal_eval(x['new_mismatch_bases'])
                x['insertion_locs'] = ast.literal_eval(x['insertion_locs'])
                x['deletion_is_micro'] = ast.literal_eval(x['deletion_is_micro'])
                x['micro_seq'] = ast.literal_eval(x['micro_seq'])
                x['count'] = int(x['count'])

            samples[summary_fn] = data

        return samples, sample_count_list


def get_mutation_event_freq(df, summary_fn):
    def assign_mutation_type(x):
        if x['num_deletions'] + x['num_insertions'] > 1:
            # More than 1 deletion or more than 1 insertion or a mix of deletions and insertions
            return 'Compound'
        elif x['num_deletions'] == 1:
            return 'Deletion'
        elif x['num_insertions'] == 1:
            return 'Insertion'
        elif x['num_mismatch']:
            return 'Mismatch'
        else:
            return 'Unclear'

    df['mutation_type'] = df.apply(assign_mutation_type, axis=1)
    count_d = {}
    count_d['Compound'] = df[df['mutation_type'] == 'Compound']['count'].sum()
    count_d['Deletion'] = df[df['mutation_type'] == 'Deletion']['count'].sum()
    count_d['Mismatch'] = df[df['mutation_type'] == 'Mismatch']['count'].sum()
    dft = pd.DataFrame(count_d.items(), columns=['mutation_type', 'count'])
    num_wt = sample_count_list[summary_fn] - dft['count'].sum()
    dft.loc[len(dft)] = ['WT', num_wt]
    dft['percentage'] = dft['count'] / dft['count'].sum()

    return dft


def del_len_dist(df):
    # get the length of all deletions
    del_lens = []
    for index, row in df.iterrows():
        for l in row['deletion_lens']:
            t_del_lens = [l] * row['count']
            del_lens.extend(t_del_lens)
    dft = pd.DataFrame({'deletion_length': del_lens})
    return dft


def seq_counts(df, mode='deletion'):

    seq_del_counts = {}
    for index, row in df.iterrows():
        count = row['count']

        if mode == 'deletion':
            seqs = re.findall('\^\(.*?\)', row['repair_sequence'])
            for seq in seqs:
                seq = seq[2:-1]
                count = seq_del_counts.get(seq, 0) + row['count']
                seq_del_counts[seq] = count

        elif mode == 'insertion':
            seqs = re.findall('\[.*?\]', row['repair_sequence'])
            for seq in seqs:
                seq = seq[1:-1]
                count = seq_del_counts.get(seq, 0) + row['count']
                seq_del_counts[seq] = count

        else:
            raise NotImplementedError

    dft = pd.DataFrame(seq_del_counts.items(), columns=['sequence', 'count'])

    return dft


def insertion_len_dist(df):
    insertion_lens = []
    for index, row in df.iterrows():
        for l in row['insertion_seqs']:
            t_ins_lens = [len(l)] * row['count']
            insertion_lens.extend(t_ins_lens)
    dft = pd.DataFrame({'insertion_length': insertion_lens})
    return dft


if __name__ == '__main__':
    # parse args
    arg_parser = argparse.ArgumentParser(description="Visualize Mutation Events")

    arg_parser.add_argument("-R", default='./sam_dsbra_plotting_ling.R',
                            help='R plotting script path')
    arg_parser.add_argument("-o", "--out_dir",
                            default='plots_%s/'%str(datetime.now())[:10],
                            help="output directory, default=summary_graph%s/"%str(datetime.now())[:10])
    arg_parser.add_argument("-i", "--input_files", nargs='+', required=True,
                            help="the output txt files from sam_dsbra_interface.py")
    arg_parser.add_argument("--all", action='store_true',
                            help="plot all types of plots")
    arg_parser.add_argument("--mut_type", action='store_true',
                            help="plot Mutation Event Frequency by Type")
    arg_parser.add_argument("--del_len", action='store_true',
                            help="plot Frequency of Deletions by Length")
    arg_parser.add_argument("--del_seq", action='store_true',
                            help="plot Sequences with Deletion Event")
    arg_parser.add_argument("--ins_len", action='store_true',
                            help="plot Frequency of Insertions by Length")
    arg_parser.add_argument("--ins_seq", action='store_true',
                            help="plot Sequences with Insertion Event")

    args = arg_parser.parse_args()

    # R script name
    r = args.R

    # read in data
    summary_fns = args.input_files

    # set output dir
    output_dir = args.out_dir
    if not os.path.exists(output_dir): os.makedirs(output_dir)

    print("Writing output to %s ..."%output_dir)

    samples, sample_count_list = parse_files(summary_fns)
    for summary_fn in summary_fns:
        summary_name = os.path.basename(summary_fn).replace('.txt', '')
        cols = [col for col in samples[summary_fn][0].keys() if col]
        df = pd.DataFrame.from_records(samples[summary_fn], columns=cols)

        # plot mutaion type frequency distribution
        if args.mut_type or args.all:
            dft = get_mutation_event_freq(df, summary_fn)
            typedist_outfn = os.path.join(output_dir,
                                          '%s_mutation_event_frequency_by_type.csv'%summary_name)
            dft.to_csv(typedist_outfn)

            # call R plot script
            typedist_outplot = os.path.join(output_dir,
                                          '%s_mutation_event_frequency_by_type.png'%summary_name)
            subprocess.call("%s --input %s --mut_type"%(r, typedist_outfn), shell=True)

        # plot deletion events: Frequency of Deletions by Length
        if args.del_len or args.all:
            dft = del_len_dist(df)
            del_freq_outfn = os.path.join(output_dir,
                                          '%s_deletion_lens.csv'%summary_name)
            dft.to_csv(del_freq_outfn)

            # call R
            del_freq_outplot = os.path.join(output_dir,
                                          '%s_deletion_lens.png'%summary_name)
            subprocess.call("%s --input %s --del_len"%(r, del_freq_outfn), shell=True)

        if args.del_seq or args.all:
            dft = seq_counts(df, mode='deletion')
            seq_with_del_outfn = os.path.join(output_dir,
                                              '%s_sequences_with_deletion_events.csv'%summary_name)
            dft.to_csv(seq_with_del_outfn)

            # call R
            seq_with_del_outplot = os.path.join(output_dir,
                                              '%s_sequences_with_deletion_events.png'%summary_name)
            subprocess.call("%s --input %s --del_seq"%(r, seq_with_del_outfn), shell=True)

        if args.ins_len or args.all:
            dft = insertion_len_dist(df)
            ins_freq_outfn = os.path.join(output_dir,
                                          '%s_insertion_lens.csv'%summary_name)
            dft.to_csv(ins_freq_outfn)

            # call R
            ins_freq_outplot = os.path.join(output_dir,
                                          '%s_insertion_lens.png'%summary_name)
            subprocess.call("%s --input %s --ins_len"%(r, ins_freq_outfn), shell=True)

        if args.ins_seq or args.all:
            dft = seq_counts(df, mode='insertion')
            seq_with_ins_outfn = os.path.join(output_dir,
                                              '%s_sequences_with_insertion_events.csv'%summary_name)
            dft.to_csv(seq_with_ins_outfn)

            # call R
            seq_with_ins_outplot = os.path.join(output_dir,
                                              '%s_sequences_with_insertion_events.png'%summary_name)
            subprocess.call("%s --input %s --ins_seq"%(r, seq_with_ins_outfn), shell=True)

        print("\nComplete!")
