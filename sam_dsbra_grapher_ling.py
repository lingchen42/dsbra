#!/usr/bin/env python
'''
This script is to visualize the mutation events in the output of sam_dsbra_interface.py
Usage:
    sam_dsbra_grapher_ling.py -h

Ling Chen 2018-05-29
'''
from __future__ import division
import os
import re
import csv
import ast
import argparse
import subprocess
import itertools
from datetime import datetime

import pandas as pd
from scipy.stats import poisson


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

    for summary_fn in summary_fns:
        with open(summary_fn,'r') as sf:

            data = list(csv.DictReader(sf, delimiter='\t'))
            #formatting of all fields appropriately
            for x in data:
                try:
                    x['insertion_seqs'] = ast.literal_eval(x['insertion_seqs'])
                    x['deletions_seqs'] = ast.literal_eval(x['deletions_seqs'])
                    x['deletion_lens'] = ast.literal_eval(x['deletion_lens'])
                    x['num_mismatch'] = int(x['num_transitions'])
                    x['num_transitions'] = int(x['num_transitions'])
                    x['num_deletions'] = int(x['num_deletions'])
                    x['num_insertions'] = int(x['num_insertions'])
                    x['new_mismatch_bases'] = ast.literal_eval(x['new_mismatch_bases'])
                    x['insertion_locs'] = ast.literal_eval(x['insertion_locs'])
                    x['deletion_is_micro'] = ast.literal_eval(x['deletion_is_micro'])
                    x['micro_seq'] = ast.literal_eval(x['micro_seq'])
                    x['repair_cigar_tuples'] = ast.literal_eval(x['repair_cigar_tuples'])
                except:  # wild type entry
                    pass

                x['repair_size'] = ast.literal_eval(x['repair_size'])
                x['count'] = int(x['count'])
                x['ref_mut_start'] = int(ast.literal_eval(x['ref_mut_start']))
                x['ref_mut_end'] = int(ast.literal_eval(x['ref_mut_end']))

            samples[summary_fn] = data

        return samples


def apply_count_cutoff(df, p=0.001, num_colony=300):
    '''
    based on poisson distribution, apply count cutoff to the data
    '''
    total_reads = df['count'].sum()  # total reads from table, not including the failed alignments.
    expected_reads = total_reads / num_colony
    count_cutoff = poisson.ppf(p, expected_reads)
    print("Using %s as count cutoff\n"%round(count_cutoff, 2))
    df = df[df['count'] >= count_cutoff]
    return df


def get_mutation_event_freq(df, summary_fn):
    def assign_mutation_type(x):
        if (x['repair_sequence']=='') or (x['repair_sequence']=='WT'):
            return "WT"
        elif x['num_deletions'] + x['num_insertions'] + x['num_mismatch'] > 1:
            # if more than 1 mutation event
            return 'Compound'
        elif x['num_deletions']:
            return 'Deletion'
        elif x['num_insertions']:
            return 'Insertion'
        elif x['num_mismatch']:
            return 'Mismatch'
        else:
            return 'Unclear'

    df['mutation_type'] = df.apply(assign_mutation_type, axis=1)
    count_d = {}
    count_d['Compound'] = df[df['mutation_type'] == 'Compound']['count'].sum()
    count_d['Deletion'] = df[df['mutation_type'] == 'Deletion']['count'].sum()
    count_d['Insertion'] = df[df['mutation_type'] == 'Insertion']['count'].sum()
    count_d['Mismatch'] = df[df['mutation_type'] == 'Mismatch']['count'].sum()
    count_d['WT'] = df[df['mutation_type'] == 'WT']['count'].sum()
    dft = pd.DataFrame(count_d.items(), columns=['mutation_type', 'count'])
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
            seqs = row['deletions_seqs']
            #seqs = re.findall('\^\(.*?\)', row['repair_sequence'])
            for seq in seqs:
                count = seq_del_counts.get(seq, 0) + row['count']
                seq_del_counts[seq] = count

        elif mode == 'insertion':
            seqs = row['insertion_seqs']
            #seqs = re.findall('\[.*?\]', row['repair_sequence'])
            for seq in seqs:
                count = seq_del_counts.get(seq, 0) + row['count']
                seq_del_counts[seq] = count

        elif mode == 'repair_seq':
            seq = row['repair_sequence']
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


def aligned_mut(df):
    df["sort"] = df['repair_size'].abs()
    df = df.sort_values("sort", ascending=True)

    indices_start = []
    indices_end = []
    mut_starts = []
    mut_ends = []
    types = []
    region_lens = []
    counts = []

    cigar_d = {1:'insertion', 2:'deletion', 7:'matched', 8:'mismatch'}
    cum_count = 0
    for index, row in df.iterrows():
        idx_start = cum_count
        idx_end = cum_count + row['count']
        cum_count += row['count']
        mut_start = row['ref_mut_start']
        mut_end = row['ref_mut_end']
        repair_cigar_tuples = row['repair_cigar_tuples']
        if repair_cigar_tuples:
            loc = mut_start
            region_len = sum(zip(*repair_cigar_tuples)[1])
            for event_type, event_len in repair_cigar_tuples:
                indices_start.append(idx_start)
                indices_end.append(idx_end)
                types.append(cigar_d[event_type])
                mut_starts.append(loc)
                counts.append(row['count'])
                region_lens.append(region_len)
                if event_type == 1:  # insertion does not consume ref bases
                    mut_ends.append(loc+0.9)  # just to make it visible
                else:
                    mut_ends.append(loc+event_len)
                    loc += event_len
        else:
            indices_start.append(idx_start)
            indices_end.append(idx_end)
            mut_starts.append(0)
            mut_ends.append(0)
            types.append("matched")
            region_lens.append(None)
            counts.append(row['count'])

    dft = pd.DataFrame()
    dft['idx_start'] = indices_start
    dft['idx_end'] = indices_end
    dft['mut_start'] = mut_starts
    dft['mut_end'] = mut_ends
    dft['type'] = types
    dft['region_len'] = region_lens
    dft['count'] = counts

    return dft


if __name__ == '__main__':
    # parse args
    arg_parser = argparse.ArgumentParser(description="Visualize Mutation Events")

    arg_parser.add_argument("-R", default='./sam_dsbra_plotting_ling.R',
                            help='R plotting script path')
    arg_parser.add_argument("-o", "--out_dir",
                            default='plots_%s/'%str(datetime.now())[:10],
                            help="output directory, default is plots_%s/"%str(datetime.now())[:10])
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
    arg_parser.add_argument("--repair_seq", action='store_true',
                            help="plot top most common repair sequences")
    arg_parser.add_argument("--aligned_mutations", action='store_true',
                            help="plot frequency and aligned mutations")
    arg_parser.add_argument("--ref_bottom", default=0, type=int,
                            help="the start reference base to show aligned mutation events")
    arg_parser.add_argument("--ref_top", default=171, type=int,
                            help="the start reference base to show aligned mutation events")
    arg_parser.add_argument("--break_index", type=int, default=None,
                            help="the index of break site")
    arg_parser.add_argument("--fts", type=int, default=14,
                            help="axis font size; default 14")
    arg_parser.add_argument("--count_cutoff", nargs=2, type=float, default=None,
                            help="take probability cutoff, number of colonies,"
                                 "apply count cutoff to the summary table"
                                 "based on the poisson distribution;"
                                 "default 0.001, 300")

    args = arg_parser.parse_args()

    # R script name
    r = args.R
    break_index = args.break_index
    ref_top = args.ref_top
    ref_bottom = args.ref_bottom

    # read in data
    summary_fns = args.input_files

    # set output dir
    output_dir = args.out_dir
    if not os.path.exists(output_dir): os.makedirs(output_dir)

    # fts
    fts = args.fts

    print("Writing output to %s ..."%output_dir)

    #samples, sample_count_list = parse_files(summary_fns)
    samples = parse_files(summary_fns)
    for summary_fn in summary_fns:
        summary_name = os.path.basename(summary_fn).replace('.txt', '')
        cols = [col for col in samples[summary_fn][0].keys() if col]
        df = pd.DataFrame.from_records(samples[summary_fn], columns=cols)
#        df = pd.read_table(summary_fn)

        # apply count cutoff
        if args.count_cutoff:
            p, num_colony = args.count_cutoff
            df = apply_count_cutoff(df, p, num_colony)

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
            subprocess.call("%s --input %s --del_seq --fts %s"%(r, seq_with_del_outfn, fts), shell=True)

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

        if args.repair_seq or args.all:
            dft = seq_counts(df, mode='repair_seq')
            repair_pattern_outfn = os.path.join(output_dir, '%s_sequences_with_mutation_events.csv'%summary_name)
            dft.to_csv(repair_pattern_outfn)
            subprocess.call("%s --input %s --repair_seq --fts %s"%(r, repair_pattern_outfn, fts), shell=True)

        if args.aligned_mutations or args.all:
            if break_index:
                dft = aligned_mut(df)
                aligned_mut_outfn = os.path.join(output_dir, '%s_aligned_mutation_events.csv'%summary_name)
                dft.to_csv(aligned_mut_outfn)
                subprocess.call("%s --input %s --aligned_mutations --break_index %s --ref_bottom %s --ref_top %s"\
                                %(r, aligned_mut_outfn, break_index, ref_bottom, ref_top), shell=True)
            else:
                print("Error generating alignment plot. Please provide break index using --break_index.")

        subprocess.call("rm %s/*csv"%output_dir, shell=True)
        print("\nComplete!")
