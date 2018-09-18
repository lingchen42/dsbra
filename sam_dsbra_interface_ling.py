'''
2018-07-31 Ling Chen
Double strand break repair analyzer (DSBRA)
Major changes
1. This script uses a customed aligner that chooses the alignment with
   mutations closest to the break site. Original aligner uses bowtie2
   and fails to get the closet alignment to the break site.
2. Output the insertion sequence and the query sequence near the repair site.
   The original script does not.
3. Also includes the count of wildtypes.

Usage: sam_dsbra_interface_ling.py [-h] -o OUT_DIR [--run_name RUN_NAME]
                                   [-r REF_FA] [-q FASTQ] [-b BREAK_INDEX]
                                   [-m MARGIN] [-qm Q_REPAIR_MARGIN]
                                   [--last_margin LAST_MARGIN]
                                   [--count_cutoff COUNT_CUTOFF COUNT_CUTOFF]
                                   [--full_table]
'''

from __future__ import division
import os
import re
import sys
import time
import shutil
import argparse
import cPickle
import tempfile
import subprocess
from itertools import groupby
from datetime import datetime
from collections import OrderedDict

import pysam
import numpy as np
import pandas as pd
from scipy.stats import poisson

SAMTOOLS_FOLDER = ''
#CIGAR_RE = re.compile('\d+[MIDNSHP]')
RE_FASTA = re.compile(r'\>(.*)[\n|\r]+([ACTGactg\r\n]+)')
TRANSITION_D = {'A':'G', 'C':'T', 'G':'A', 'T':'C'}

class RepairSample:
    def __init__(self):
        self.total_samples

class RepairPattern:
    def __init__(self, read, read_name, reference_end):
        self.read = read
        self.read_name = read_name
        self.repair_size = 0

        self.repair_sequence = ''
        self.repair_cigar_tuples = []
        # repair seq neighbourhood (+/- q_repair_margin on either side of the  mut range)in query seq
        self.actual_seq_near_break = ''
        self.mut_bottom_range = None
        self.mut_top_range = None
        self.ref_mut_start = 0
        self.ref_mut_end  = reference_end

        self.mismatch_locations = [ ]
        self.new_mismatch_bases = [ ]
        self.old_mismatch_bases = [ ]
        self.num_transitions = 0
        self.num_transversions = 0

        self.num_insertions = 0
        self.insertion_locs = [ ]
        self.insertion_seqs = [ ]

        self.num_deletions = 0
        self.deletion_locs = [ ]
        self.deletion_lens = [ ]
        self.deletions_seqs = [ ]
        self.deletion_is_micro = [ ]
        self.micro_seq = [ ]

    def add_insertion(self, ins_loc, ins_str):
        self.num_insertions += 1
        self.insertion_locs.append(ins_loc)
        self.insertion_seqs.append(ins_str)

    def add_deletion(self, del_loc, del_size, del_str, is_micro, micro_seq):
        self.num_deletions += 1
        self.deletion_locs.append(del_loc)
        self.deletion_lens.append(del_size)
        self.deletions_seqs.append(del_str)
        self.deletion_is_micro.append(is_micro)
        self.micro_seq.append(micro_seq)

    def add_mismatch(self, mismatch_loc, new_base, old_base):
        self.mismatch_locations.append(mismatch_loc)
        self.new_mismatch_bases.append(new_base)
        self.old_mismatch_bases.append(old_base)

        if TRANSITION_D[old_base] == new_base:
            self.num_transitions += 1
        else:
            self.num_transitions += 1

    def get_repair_sequence(self, read, ref_seq, q_repair_margin):
        '''
        Generate repair sequences
        deletion: ^()
        insertion: []
        mismatch: mismatched base A/G/C/T
        matched: *
        '''
        repair_seq = self.repair_sequence

        q_seq = read.query_sequence
        q_aln_locs, r_aln_locs = zip(*read.get_aligned_pairs())  # the length of these two should be equal to the sum of cigar
        q_aln_locs = list(q_aln_locs)
        r_aln_locs = list(r_aln_locs)
        mut_start = self.mut_bottom_range
        mut_end = self.mut_top_range

        # mut_start/end in ref for repair seq
        self.ref_mut_start = r_aln_locs[mut_start]
        if self.ref_mut_start == None:  # in case of insertion, it will be None; don't use if not, because 0 is the start index
            self.ref_mut_start = r_aln_locs[mut_start-1]

        if mut_end < len(r_aln_locs):  # incase the deletion extends to the end
            self.ref_mut_end = r_aln_locs[mut_end]
        else:
            self.ref_mut_end = r_aln_locs[mut_end-1]

        # mut_start/end in reads
        q_repair_start = q_aln_locs[mut_start]
        if q_aln_locs[mut_start] == None:  # in case of deletion it will be NaN
            if (mut_start-1) > 0:
                q_repair_start = q_aln_locs[mut_start-1]
            else:
                q_repair_start = 0

        if mut_end >= len(ref_seq):  # in case the deletion extends to the end
            q_repair_end = q_aln_locs[len(ref_seq) - 1]
        else:
            q_repair_end = q_aln_locs[mut_end]
        if  q_repair_end == None:  # in case of deletion it will be NaN
            if ((mut_end+1) < len(q_aln_locs)):
                q_repair_end = q_aln_locs[mut_end+1]
            else:
                # in case the deletion extends to the end
                q_repair_end = len(q_seq) - 1

        q_repair_start_m = max(q_repair_start - q_repair_margin, 0)
        q_repair_end_m = min(q_repair_end + q_repair_margin, len(q_seq))
        self.actual_seq_near_break = q_seq[q_repair_start_m : q_repair_end_m]

        # cigar_d = {1:"insertion", 2:"deletion", 7:"matched", 8:"mismatched"}  # cannot deal with other status
        cigar_tuples = read.cigartuples
        cigar_list = []
        for event_type, event_len in cigar_tuples:
            cigar_list.extend([event_type for i in range(event_len)])
        cigar_list = cigar_list[mut_start : mut_end]
        groups = groupby(cigar_list)
        cigar_tuples_in_region = [(label, sum(1 for _ in group)) for label, group in groups]

        loc = mut_start
        for event_type, event_len in cigar_tuples_in_region:
            q_loc = q_aln_locs[loc]
            r_loc = r_aln_locs[loc]
            if event_type == 1:  # insertion
                repair_seq += "[%s]"%q_seq[q_loc : q_loc + event_len]
                self.repair_size += event_len
            elif event_type == 2:  # deletion
                repair_seq += "^(%s)"%ref_seq[r_loc : r_loc + event_len]
                self.repair_size -= event_len
            elif event_type == 7:  # match
                repair_seq += "*"*event_len
            elif event_type == 8:  # mismatch
                repair_seq += q_seq[q_loc : q_loc + event_len]
            else:
                print("Error! Unexpected cigar number %s detected"%event_type)
                sys.exit(1)
            loc += event_len

        self.repair_cigar_tuples = cigar_tuples_in_region
        self.repair_sequence = repair_seq


def get_ref(ref_fa, output_dir, run_name):
    with open(ref_fa, 'r') as ref_fasta:
        ref_fasta = open(ref_fa,'r')
        ref_fasta_str = ref_fasta.read()
        index_ref_fn = os.path.join(output_dir, '%s_ref.fa'%run_name)
        with open(index_ref_fn, 'w') as index_ref:
            first_seq = RE_FASTA.match(ref_fasta_str)
            ref_name = ''
            ref_seq = ''
            if first_seq:
                ref_name, ref_seq = first_seq.groups(0)
                ref_seq = re.sub(r'[\n\r]','',ref_seq).upper()
                index_ref.write('>'+ref_name+'\n'+ref_seq+'\n')
            else:
                raise EOFError('FASTA sequence not detected before end of file.'\
                               + str(ref_fa))
    return ref_seq


def aln(sam_filename, ref_fa, fastq_name, break_index,
        pcr_primer1, pcr_primer2, check_n_nuc_at_end,
        pcr_tail_seq, min_len, not_valid_reads,
        output_dir, df_run_info,
        run_info, align_mode, record_failed_alignment=True):
    '''
    Customed aligner. Among alignments with highest scores, choose the one with
    mutation events cloest to the break site.
    Args:
        sam_filename: the output alignment filename
        ref_fa: the path to the reference fasta file
        fastq_name: the path to the read fastq file
        break_index: where the break occurs. O-indexed, the position of the
                     nucleotide right before the break site.
        df_run_info: For recording alignment summary
        record_failed_alignment: Whether write the failed alignment to a file.
    Returns:
        run_metadata: run_metadata dictionary with alignment summary
        It will write the alignments to sam_filename.
    '''

    # perform the alignment, sorting, indexing
    with open(os.devnull, 'wb') as out:
        if not os.path.exists(sam_filename):  # if no corresponding SAM file, do the alignment
            print("\nStart alignment...")
            aln_cmd = './custom_alignment_dsbra.py -o %s '\
                      '--run_info %s --ref %s '\
                      '-f %s -b %s '\
                      '--align_mode %s '\
                      '--pcr_primer1 %s --pcr_primer2 %s '\
                      '--pcr_tail_seq %s --min_len %s '\
                      %(sam_filename, run_info, ref_fa, fastq_name, break_index,
                        align_mode, pcr_primer1, pcr_primer2,
                        pcr_tail_seq, min_len)
            if not_valid_reads:
                aln_cmd = aln_cmd + " --not_valid_reads %s"%not_valid_reads
            df_run_info['alignment_settings'] = aln_cmd
            # save the run info because custom alignment will need to read it
            df_run_info.to_csv(run_info, index=False)
            print(df_run_info['alignment_settings'][0])
            subprocess.call(df_run_info['alignment_settings'], shell=True)
            print("Alignment completed. The alignment file is:\n%s\n"%sam_filename)
            # read back the dataframe
            df_run_info = pd.read_csv(run_info)
        else:
            print("Alignment file %s already exist. Will use the existing alignment file."%sam_filename)

        subprocess.call(SAMTOOLS_FOLDER+'samtools view -@ 8 -bS '+sam_filename+' | '+SAMTOOLS_FOLDER+'samtools sort -@ 8 -o '\
                        +sam_filename+'.sorted.bam', shell=True, stdout = out, stderr=subprocess.STDOUT)
        subprocess.call(SAMTOOLS_FOLDER+'samtools index '+ sam_filename + '.sorted.bam ' + sam_filename + '.sorted.bai',\
                        shell=True, stdout = out)

    # output the failed alignment to a file
    if record_failed_alignment:
        failed_alignment_fn = os.path.join(output_dir,
                     sam_filename.split("/")[-1][:-4] + '_failed_alignment.sam')
        print("Writing failed_alignment to %s"%failed_alignment_fn)
        failed_aln = subprocess.call("awk 'FNR > 3 {if($2==4){print $0}}' "\
                                    + sam_filename + '>' +  failed_alignment_fn,
                                     shell=True)

    # summary of alignment
    cmd = SAMTOOLS_FOLDER +'samtools flagstat %s'%(sam_filename + '.sorted.bam')
    aln_stats = subprocess.check_output(cmd, shell=True)
    n_total = int(aln_stats.split('\n')[0].split(' + ')[0])
    n_success_aln = int(aln_stats.split('\n')[4].split(' + ')[0])
    n_failed_aln = n_total - n_success_aln
    df_run_info['num_total_aligned_reads'] = n_total
    df_run_info['num_failed_alignments'] = n_failed_aln
    df_run_info['num_success_alignments'] = n_success_aln
    print("success alignments: %s"%n_success_aln)

    return df_run_info


def check_microhomo_one_way(ref_seq, del_start_loc, del_len, m, mode='fw'):
    is_mmej = False
    microhomology = ''
    matched = True
    while matched:
        if mode == 'fw':  # check fw microhomology
            inside_del_span = (del_start_loc, del_start_loc + m)
            outside_del_span = (del_start_loc + del_len, del_start_loc + del_len + m)
        else:  # check reverse microhomology
            inside_del_span = (del_start_loc + del_len - m, del_start_loc + del_len)
            outside_del_span = (del_start_loc - m, del_start_loc)

        if ref_seq[inside_del_span[0] : inside_del_span[1]] \
           == ref_seq[outside_del_span[0] : outside_del_span[1]]:
            m += 1
            is_mmej = True
            microhomology = ref_seq[inside_del_span[0] : inside_del_span[1]]
        else:
            matched = False
    return is_mmej, microhomology



def check_microhomology(ref_seq, del_start_loc, del_len, microhomology_cutoff=1):
    '''
    check is the deletion is microhomology
    Args:
        ref_seq: reference sequence
        del_end_loc: the position of the last nucleotide of the deletion in reference sequence
        del_len: length of deletion
        microhomology_cutoff: the minimum length to count as being microhomology
    Return:
        is_mmej: boolean, whether it's a microhomology or not
        microhomology: microhomology sequence
    '''
    m = microhomology_cutoff
    is_mmej, microhomology = check_microhomo_one_way(ref_seq, del_start_loc,
                                                     del_len, m, mode='fw')
    if not is_mmej:   # if not forward microhomology, chech reverse.
        is_mmej, microhomology = check_microhomo_one_way(ref_seq, del_start_loc,
                                                         del_len, m, mode='rc')
    return is_mmej, microhomology


def extend_for_compound(read, rp, last_margin, mut_range):
    mut_start = min(mut_range)
    mut_end = max(mut_range)
    # for compound repair patterns
    if rp.num_deletions + rp.num_insertions + rp.num_transitions \
       + rp.num_transversions > 1:
        cigar_list = []
        for event_type, event_len in read.cigartuples:
            cigar_list.extend([event_type]*event_len)

        d = 1
        while d <= last_margin and ((mut_start-d) >= 0):
            if cigar_list[mut_start - d] != 7:  # not matched
                mut_start -= d  # update the mut_start
                d = 1  # reset distance to the last mutation events
            else:
                d += 1  # look 1 bp far away from the mut_start
        d = 1
        while d <= last_margin and ((mut_end+d) < len(cigar_list)):
            if cigar_list[mut_end + d] != 7:  # not matched
                mut_end += d  # update the mut_end
                d = 1  # reset distance to the last mutation events
            else:
                d += 1  # look 1 bp far away from the mut_end

    return mut_start, mut_end


def mutation_pattern(read, ref_seq, break_index, margin, last_margin,
                     mismatch_cutoff=1):
    '''
    create a repair pattern entry for the read
    Args:
        read: pysam.AlignedSegment object
        ref_seq: reference DNA sequence
        break_index:
        margin:
        last_margin:
        mismatch_cutoff: the mismatch event has to be < N bp from the break site
                         to be considered in the repair size/pattern calculation
    Returns:
        rp: RepairPattern class object
    '''

    rp = RepairPattern(read, read.query_name, read.reference_end)

    q_aln_locs, r_aln_locs = zip(*read.get_aligned_pairs())  # the length of these two should be equal to the sum of cigar
    q_aln_locs = list(q_aln_locs)
    r_aln_locs = list(r_aln_locs)

    # map the region of interest based on the break index and margin
    global_locs = range(len(r_aln_locs))
    region_start = r_aln_locs.index(break_index + 1 - margin)  # if the break index is 84, we want [75, 95)
    region_end = r_aln_locs.index(break_index + 1 + margin)

    # if the mutation event overlaps with region of interest, record it.
    cigar_tuples = read.cigartuples
    loc = 0
    mut_range = []  # based on global locs
    for event_type, event_len in cigar_tuples:
        q_loc = q_aln_locs[loc]  # map loc to location in query sequence
        r_loc = r_aln_locs[loc]  # map loc to location in reference sequence

        # if there is any overlap of mutation event with the region of interest,
        # we'll analyze it
        overlap = set(range(loc, loc + event_len)) & set(range(region_start, region_end))
        if overlap:
            if event_type == 7:  # matched
                pass

            elif event_type == 8:  # mismatch
                # only consider mismatches that is really close to the breaksite
                if abs(break_index - r_loc) <= mismatch_cutoff:
                    mut_range.extend([loc, loc + event_len])
                    for l in range(loc, loc + event_len):
                        rp.add_mismatch(l, read.query_sequence[q_aln_locs[l]],
                                        ref_seq[r_aln_locs[l]])

            elif event_type == 1:    # insertion
                mut_range.extend([loc, loc + event_len])
                ins_str = read.query_sequence[q_loc : q_loc + event_len]
                if loc > 0:  # if the insertion does not occur before ref starts
                    rp.add_insertion(r_aln_locs[loc - 1], ins_str)

            elif event_type == 2:  # deletion
                mut_range.extend([loc, loc + event_len])
                del_str = ref_seq[r_loc : r_loc + event_len]
                is_mmej, microhomology = check_microhomology(ref_seq, r_loc, event_len)
                rp.add_deletion(r_loc, event_len, del_str, is_mmej, microhomology)

            else:
                print("Cigar type %s is not implemented. Skip this read: %s"\
                      %(event_type, cigar_tuples))
                return "skipped"

        loc += event_len  # update current loc

    # if it's a compound mutation (not a single deletion or insertion),
    # then use the last margin to extend the region to find repair pattern
    # in a broader range. This is good for detecting telomere-like seq.
    if mut_range:
        rp.mut_bottom_range, rp.mut_top_range = extend_for_compound(read, rp, last_margin, mut_range)
    else:
    # if the mismatch > mismatch_cutoff, but there are no other mutation events within break +/- margin
    # then the mut_bottom_range or mut_top_range may have not been changed. However, it's likely
    # sequencing error, we are going to ignore this read.
        return ("WT", rp)

    return ("Mut", rp)


def get_rp(read, ref_seq, rp_entries, break_index, margin, last_margin,
           q_repair_margin):
        result = mutation_pattern(read, ref_seq, break_index, margin, last_margin)
        rp = result[1]
        if result[0] != 'WT':
            rp.get_repair_sequence(read, ref_seq, q_repair_margin)
        else:
            rp.repair_sequence = 'WT'

        for key in rp_entries.keys():
             if key != "num_mismatch":
                 rp_entries[key].append(getattr(rp, key))
        rp_entries['num_mismatch'].append((rp.num_transitions\
                                          + rp.num_transversions))

        return rp_entries


def get_summary_table(repair_pattern_table_fn, df_run_info):

    # Reading the table from file solves the problem caused by a python list in
    # a cell. Otherwise it cannot be grouped appropriately
    df = pd.read_table(repair_pattern_table_fn)
    df = df.drop("read_name", axis=1)
    df['count'] = 1
    cols_order = df.columns
    # In case actual_seq_near_break are different when repair_sequence is the same
    cols_to_grp = [c for c in cols_order if c not in ['count', 'actual_seq_near_break']]
    df = df.groupby(cols_to_grp, as_index=False).agg({'count':'sum',
                                                      'actual_seq_near_break':'first'})
    df['mut_freq'] = df['count'] / (df['count'].sum())
    df = df.sort_values(by='count', ascending=False)
    df = df[cols_order]

    return df


def apply_count_cutoff(df, df_run_info, p=0.001, num_colony=300):
    '''
    based on poisson distribution, apply count cutoff to the data
    '''
    total_reads = int(df_run_info['num_failed_alignments'][0]) \
                  + int(df_run_info['num_success_alignments'][0])
    expected_reads = total_reads / num_colony
    count_cutoff = poisson.ppf(p, expected_reads)
    print("Using %s as count cutoff\n"%round(count_cutoff, 2))
    df = df[df['count'] >= count_cutoff]
    return df


def main():
    # parse args
    arg_parser = argparse.ArgumentParser(description="Double strand break repair analyzer (DSBRA)")

    arg_parser.add_argument("-o", "--out_dir", required=True,
                            help="output directory")
    arg_parser.add_argument("--align_mode", default="global", choices=["global",
                            "local"], help="choose the method for alignment")
    arg_parser.add_argument("--pcr_primer1", required=True,
                            help="primer seq, at the begining of the ref seq; will require the read to starts with this primer seq")
    arg_parser.add_argument("--pcr_primer2", required=True,
                            help="primer seq, at the end of the ref seq; will require last 4 nucleotide of the read is part of this primer seq")
    arg_parser.add_argument("--pcr_tail_seq", required=True,
                            help="added sequence for PCR;default ATCGGAAGAGCACACGTCTGAACTCCAGTCAC")
    arg_parser.add_argument("--min_len", type=int, default=44,
                            help="miminal valid length (reads - pcr_tail_seq); default 44bp")
    arg_parser.add_argument("--not_valid_reads", default=False,
                            action='store_true', help="output not valid reads")
    arg_parser.add_argument("-r", "--ref_fa",
                            help="path of the reference sequence fasta file")
    arg_parser.add_argument("-q", "--fastq",
                            help="path of the reads to analyze in fastq format")
    arg_parser.add_argument("-b", "--break_index", type=int,
                            help="break index. 0-indexed, the position of"
                                 "nucleotide right before the break site.")
    arg_parser.add_argument("-m", "--margin", type=int, default=5,
                            help="the margin before and after the break index"
                                 "that should be considered for analysis;"
                                 "default 5bp.")
    arg_parser.add_argument("-qm", "--q_repair_margin", type=int, default=5,
                            help="the margin before and after the last mutation"
                                 "that should be reported"
                                 "in the repaired_read_sequence; default 5bp.")
    arg_parser.add_argument("--last_margin", type=int, default=5,
                            help="the margin outside of the last mutation event"
                                "that should be included in the repair sequence"
                                "pattern for those compound mutation events"
                                "; default 5bp")
    arg_parser.add_argument("--count_cutoff", nargs=2, type=float, default=None,
                            help="take probability cutoff, number of colonies,"
                                 "apply count cutoff to the summary table"
                                 "based on the poisson distribution;"
                                 "default 0.001, 300")
    arg_parser.add_argument("--full_table", action="store_true", default=True,
                            help="Will generate a large table with"
                                 "repair pattern for each analyzed")
#    arg_parser.add_argument( "--intermediate_files", action="store_true",
#                            help="Whether to keep the intermediate files generated by the script")

    print("Start..")
    args = arg_parser.parse_args()
    start_time = time.time()
    print("")

    # set global values
    output_dir = args.out_dir
    if os.path.exists(output_dir) :
        print("WARNING! Overwritting the existing output directory %s."%(output_dir))
    else:
        os.makedirs(output_dir)

    align_mode = args.align_mode
    break_index = args.break_index
    margin = args.margin
    fastq_name = args.fastq
    ref_fa = args.ref_fa
    run_name = fastq_name.split("/")[-1].replace('.gz', '').split(".")[0] \
               + '_b%s_m%s'%(break_index, margin)
    sam_filename = os.path.join(output_dir, run_name + '_%s_aln.sam'%align_mode)
    q_repair_margin = args.q_repair_margin
    last_margin = args.last_margin

    summary_table_fn = os.path.join(output_dir, "%s_repair_pattern_summary_table.txt"%(run_name))
    if args.full_table:
        repair_pattern_table_fn = os.path.join(output_dir, \
                                               "%s_repair_pattern_full_table.txt"%(run_name))

    output_filename = os.path.join(output_dir, "%s_repair_pattern_summary.txt"%run_name)
    run_info = os.path.join(output_dir, '%s_run_info.csv'%run_name)

    # read in reference fasta, write to new file after checking filetype
    ref_seq = get_ref(ref_fa, output_dir, run_name)

    # record run info
    try:
        df_run_info = pd.read_csv(run_info)
    except:
        df_run_info = pd.DataFrame()
    df_run_info['fastq_file'] = fastq_name
    df_run_info["time"] = [str(datetime.now())[:-7]]
    df_run_info["working_dir"] = os.getcwd()
    df_run_info["commad"] = ' '.join(sys.argv)
    df_run_info['break_index'] = break_index
    df_run_info['margin'] = margin
    df_run_info["ref_filename"] = ref_fa
    df_run_info['ref_seq'] = ref_seq

    # customed aligner
    # whether to print not valid reads
    not_valid_reads = args.not_valid_reads
    if not_valid_reads:
        not_valid_reads =  os.path.join(output_dir,
                                        run_name + '_%s_not_valid_reads.txt'\
                                        %align_mode)
    # params for removing pcr tail seq
    pcr_tail_seq = args.pcr_tail_seq
    min_len = args.min_len
    # params for checking primer existance
    pcr_primer1 = args.pcr_primer1
    pcr_primer2 = args.pcr_primer2
    check_n_nuc_at_end = 4
    df_run_info = aln(sam_filename, ref_fa, fastq_name, break_index,
                      pcr_primer1, pcr_primer2, check_n_nuc_at_end,
                      pcr_tail_seq, min_len, not_valid_reads,
                      output_dir, df_run_info, run_info, align_mode,
                      record_failed_alignment=True)
    print("")

    # analyzing the repair pattern of each read, store in a table
    print("Analyzing repair pattern...")
    bamfile = pysam.AlignmentFile(sam_filename+".sorted.bam","rb")

    fields = ['read_name', 'repair_size', 'mismatch_locations',
              'new_mismatch_bases', 'num_mismatch', 'num_transitions',
              'num_transversions', 'num_insertions', 'insertion_locs',
              'insertion_seqs', 'num_deletions', 'deletion_locs',
              'deletion_lens', 'deletions_seqs',
              'deletion_is_micro', 'micro_seq', 'repair_sequence',
              'actual_seq_near_break', 'repair_cigar_tuples',
              'ref_mut_start', 'ref_mut_end']
    rp_entries = {}
    for f in fields:
        rp_entries[f] = []

    for read in bamfile.fetch(bamfile.references[0]):
        rp_entries = get_rp(read, ref_seq, rp_entries, break_index, margin,
                            last_margin, q_repair_margin)

    df = pd.DataFrame(rp_entries, columns=fields)
    if args.full_table:
        print("Writing full table to %s"%repair_pattern_table_fn)
        df.to_csv(repair_pattern_table_fn, sep="\t", index=False)

    # generate summary table of reads. Groups reads with the same repair pattern
    df = get_summary_table(repair_pattern_table_fn, df_run_info)
    # the reads made to the end (repair pattern analyzed reads) can be different
    # from reads_with_mutations_within_range
    # because some reads can only mismatches within range, but not immediately close to break index and got filtered out
    df_run_info['Repair pattern analyzed reads'] = df['count'].sum()
    if df[df['repair_sequence']=='WT']['count'].values:
        df_run_info['Number of wild type'] = int(df[df['repair_sequence']=='WT']['count'])
    else:
        df_run_info['Number of wild type'] = 0
    df_run_info['Number of mutated reads'] =\
        df_run_info['Repair pattern analyzed reads'][0]\
        - df_run_info['Number of wild type'][0]

    # apply count cutoff
    if args.count_cutoff:
        print("\nApply count cutoff...")
        df = apply_count_cutoff(df, df_run_info,
                                p=args.count_cutoff[0],
                                num_colony=args.count_cutoff[1])

    print("Writing summary table to %s"%summary_table_fn)
    df.to_csv(summary_table_fn, sep='\t', index=False)
    summary_table_excel = summary_table_fn.replace(".txt", ".xlsx")
    print("Writing summary table in excel format to %s \n"%summary_table_excel)
    df = df[list(df.columns[:-4]) + ['count']]
    writer = pd.ExcelWriter(summary_table_excel)
    df.to_excel(writer,'Sheet1', index=False)
    writer.save()

    # write run metadata
    print("Writing run metadata to %s\n"%(run_info))
    df_run_info.to_csv(run_info, index=False)

    print("Complete!")
    print "Total time: " + str(time.time()-start_time)


if __name__ == '__main__':
    main()
