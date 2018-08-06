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
    def __init__(self, read_name, reference_end):
        self.read_name = read_name
        self.repair_size = 0

        self.repair_sequence = ''
        self.repair_cigar_tuples = []
        # repair seq neighbourhood (+/- q_repair_margin on either side of the  mut range)in query seq
        self.actual_seq_near_break = ''
        self.mut_bottom_range = reference_end  # for the first comparison
        self.mut_top_range = -1  # for the first comparison
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

    def get_repair_sequence(self, read, ref_seq):
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
        if not self.ref_mut_start:  # in case of insertion, it will be NaN
            self.ref_mut_start = r_aln_locs[mut_start-1]
        self.ref_mut_end = r_aln_locs[mut_end]

        # mut_start/end in reads
        q_repair_start = q_aln_locs[mut_start]
        if not q_aln_locs[mut_start]:  # in case of deletion it will be NaN
            q_repair_start = q_aln_locs[mut_start-1]
        q_repair_end = q_aln_locs[mut_end]
        if not q_aln_locs[mut_end]:  # in case of deletion it will be NaN
            q_repair_end = q_aln_locs[mut_end+1]
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
                self.repair_size += event_len
            elif event_type == 8:  # mismatch
                repair_seq += q_seq[q_loc : q_loc + event_len]
                self.repair_size += event_len
            else:
                print("Error! Unexpected cigar number %s detected"%event_type)
                sys.exit(1)
            loc += event_len

        self.repair_cigar_tuples = cigar_tuples_in_region
        self.repair_sequence = repair_seq


def get_ref(ref_fa):
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


def aln(sam_filename, ref_fa, fastq_name, break_index, run_metadata, output_dir,
        record_failed_alignment=True):
    '''
    Customed aligner. Among alignments with highest scores, choose the one with
    mutation events cloest to the break site.
    Args:
        sam_filename: the output alignment filename
        ref_fa: the path to the reference fasta file
        fastq_name: the path to the read fastq file
        break_index: where the break occurs. O-indexed, the position of the
                     nucleotide right before the break site.
        run_metadata: run_metadata dictionary. For recording alignment summary
        record_failed_alignment: Whether write the failed alignment to a file.
    Returns:
        run_metadata: run_metadata dictionary with alignment summary
        It will write the alignments to sam_filename.
    '''

    # perform the alignment, sorting, indexing
    with open(os.devnull, 'wb') as out:
        if not os.path.exists(sam_filename):  # if no corresponding SAM file, do the alignment
            print("\nStart alignment...")
            run_metadata['alignment_settings'] = './custom_alignment_dsbra.py -o %s --ref %s -f %s -b %s'%(sam_filename, ref_fa, fastq_name, break_index)
            print(run_metadata['alignment_settings'])
            subprocess.call(run_metadata['alignment_settings'], shell=True)
            print("Alignment completed. The alignment file is:\n%s\n"%sam_filename)
        else:
            print("Alignment file %s already exist. Will use the existing alignment file."%sam_filename)

        subprocess.call(SAMTOOLS_FOLDER+'samtools view -@ 8 -bS '+sam_filename+' | '+SAMTOOLS_FOLDER+'samtools sort -@ 8 -o '\
                        +sam_filename+'.sorted.bam', shell=True, stdout = out, stderr=subprocess.STDOUT)
        subprocess.call(SAMTOOLS_FOLDER+'samtools index '+ sam_filename + '.sorted.bam ' + sam_filename + '.sorted.bai',\
                        shell=True, stdout = out)

    # output the failed alignment to a file
    if record_failed_alignment:
        failed_alignment_fn = os.path.join(output_dir, sam_filename.split("/")[-1][:-4] + '_failed_alignment.sam')
        print("Writing failed_alignment to %s"%failed_alignment_fn)
        failed_aln = subprocess.call("awk 'FNR > 3 {if($2==4){print $0}}' "\
                                     + sam_filename + '>' +  failed_alignment_fn,
                                     shell=True)

    # summary of alignment
    n_failed_aln = subprocess.check_output('''awk 'FNR > 3 {if($2==4){failed += 1}} END{print failed}' '''\
                                          + sam_filename, shell=True)
    n_success_aln = subprocess.check_output('''awk 'FNR > 3 {if($2!=4){s += 1}} END{print s}' '''\
                                          + sam_filename, shell=True)
    run_metadata['num_failed_alignments'] = n_failed_aln.strip()
    run_metadata['num_success_alignments'] =  n_success_aln.strip()

    return run_metadata


def mut_within_margin(sam_filename, mut_bamfile_name, non_mut_samfile_name,
                      break_index, margin, run_metadata):
    '''
    Reads with mutation within margin.
    Args:
        sam_filename: alignment file
        mut_bamfile_name: mutated reads bamfile name
        non_mut_samfile_name: not mutated reads bamfile name
        run_metadata: run_metadata dictionary. For recording average coverage of
                      bases within margin and number of qualified reads.
    Returns:
        run_metadata
        It will also write the qualified reads to a file.

    '''

    bamfile = pysam.AlignmentFile(sam_filename+'.sorted.bam','rb')

    # VERSION 0.0.4 -- Mutational database now incorporates
    # First pass on SAM/BAM file creates secondary file
    # with all reads requiring further evaluation (non-WT within margin)
    reads_list = [ ]
    mut_reads_list = [ ]

    if os.path.exists(mut_bamfile_name):
        print "ERROR: File \"%s\" exists."%mut_bamfile_name
        sys.exit(1)

    mut_bamfile = pysam.AlignmentFile(mut_bamfile_name,"wh", template=bamfile,
                                      reference_names=bamfile.references,
                                      reference_lengths=bamfile.lengths)

    print("...\nBeginning first scan of SAM file...")
    start_time = time.time()

    #Find average coverage within margin
    coverage_list = []

    for pileupcolumn in bamfile.pileup(bamfile.references[0],
                                       start=break_index-margin+1, # break index 84, we want [75, 95)
                                       end=break_index+margin+1,
                                       truncate=True,
                                       max_depth=250000):
        coverage_list.append(pileupcolumn.n)

        for pileupread in pileupcolumn.pileups:
            qname = pileupread.alignment.query_name
            reads_list.append(qname)

            if qname in mut_reads_list:
                continue

            mutated = False

            # is there is indel
            if pileupread.indel != 0 or pileupread.is_del:
                mutated = True
            # is there is a mismatch
            else:
                if ref_seq[pileupcolumn.pos] != \
                   pileupread.alignment.query_sequence[pileupread.query_position]:
                    mutated = True

            if mutated:
                mut_reads_list.append(qname)
                mut_bamfile.write(pileupread.alignment)

    non_mut_reads_list = list(set(reads_list) - set(mut_reads_list))
    non_mut_fh, non_mut_path = tempfile.mkstemp()
    try:
        with os.fdopen(non_mut_fh, 'w') as tmp:
            print("Writing reads without mutations within range to file %s"%non_mut_samfile_name)
            for r in non_mut_reads_list:
                tmp.write(r+"\n")
        cmd = 'grep -f %s %s > %s'%(non_mut_path, sam_filename, non_mut_samfile_name)
        subprocess.call(cmd, shell=True)
    finally:
        os.remove(non_mut_path)

    run_metadata['coverage'] = int(np.mean(coverage_list))
    run_metadata['num_reads_with_mutations_within_range'] = len(mut_reads_list)
    run_metadata['num_reads_without_mutations_within_range'] = len(non_mut_reads_list)

    #Close and sort/index new bamfile
    mut_bamfile.close()
    with open(os.devnull, 'wb') as out:
        subprocess.call(SAMTOOLS_FOLDER+'samtools view -@ 8 -bS '+mut_bamfile_name+' | '+SAMTOOLS_FOLDER+'samtools sort -@ 8 -o '\
            +mut_bamfile_name+'.sorted.bam', shell=True, stdout = out)
        subprocess.call(SAMTOOLS_FOLDER+'samtools index '+ mut_bamfile_name + '.sorted.bam ' + mut_bamfile_name + '.sorted.bai',\
            shell=True, stdout = out)
    bamfile.close()

    print "First scan and sorting ended. Total time: " + str(time.time()-start_time) + "\n..."

    return run_metadata


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



def check_microhomology(ref_seq, del_start_loc, del_len, microhomology_cutoff=2):
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
        while d <= last_margin:
            if cigar_list[mut_start - d] != 7:  # not matched
                mut_start -= d  # update the mut_start
                d = 1  # reset distance to the last mutation events
            else:
                d += 1  # look 1 bp far away from the mut_start
        d = 1
        while d <= last_margin:
            if cigar_list[mut_end + d] != 7:  # not matched
                mut_end += d  # update the mut_end
                d = 1  # reset distance to the last mutation events
            else:
                d += 1  # look 1 bp far away from the mut_end

    return mut_start, mut_end


def mutation_pattern(read, ref_seq, break_index, margin, mismatch_cutoff=1):
    '''
    create a repair pattern entry for the read
    Args:
        read: pysam.AlignedSegment object
        ref_seq: reference DNA sequence
        break_index:
        margin:
        mismatch_cutoff: the mismatch event has to be < N bp from the break site
                         to be considered in the repair size/pattern calculation
    Returns:
        rp: RepairPattern class object
    '''

    rp = RepairPattern(read.query_name, read.reference_end)

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
                        rp.add_mismatch(loc, read.query_sequence[q_aln_locs[l]],
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

    mut_range = [i for i in mut_range if i]

    # if it's a compound mutation (not a single deletion or insertion),
    # then use the last margin to extend the region to find repair pattern
    # in a broader range. This is good for detecting telomere-like seq.
    if mut_range:
        rp.mut_bottom_range, rp.mut_top_range = extend_for_compound(read, rp, last_margin, mut_range)
    else:
    # if the mismatch > mismatch_cutoff, but there are no other mutation events within break +/- margin
    # then the mut_bottom_range or mut_top_range may have not been changed. However, it's likely
    # sequencing error, we are going to ignore this read.
       return "skipped"

    return rp


def get_rp(read, ref_seq, rp_entries, break_index, margin):
        rp = mutation_pattern(read, ref_seq, break_index, margin)
        if rp != 'skipped':
            rp.get_repair_sequence(read, ref_seq)
            for key in rp_entries.keys():
                if key != "num_mismatch":
                    rp_entries[key].append(getattr(rp, key))
            rp_entries['num_mismatch'].append((rp.num_transitions\
                                               + rp.num_transversions))

        return rp_entries


def get_summary_table(repair_pattern_table_fn, run_metadata,
                      non_mut_samfile_name):

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

    # add entry for WT
    wt_entry = pd.Series()
    for f in list(df.columns[:-1]):
        wt_entry[f] = ""
    wt_entry['repair_size'] = 0
    wt_entry['ref_mut_start'] = 0
    wt_entry['ref_mut_end'] = 0
    try:
        wt_entry['count'] = run_metadata["num_reads_without_mutations_within_range"]
    except KeyError:
        cmd = "wc -l < %s"%non_mut_samfile_name
        wt_count = int(subprocess.check_output(cmd, shell=True))
        wt_entry['count'] = wt_count
    df = df.append(wt_entry, ignore_index=True)

    df['mut_freq'] = df['count'] / (df['count'].sum())
    df = df.sort_values(by='count', ascending=False)
    df = df[cols_order]

    return df


def apply_count_cutoff(df, run_metadata, p=0.001, num_colony=300):
    '''
    based on poisson distribution, apply count cutoff to the data
    '''
    total_reads = int(run_metadata['num_failed_alignments']) \
                  + int(run_metadata['num_success_alignments'])
    expected_reads = total_reads / num_colony
    count_cutoff = poisson.ppf(p, expected_reads)
    print("Using %s as count cutoff\n"%round(count_cutoff, 2))
    df = df[df['count'] >= count_cutoff]
    return df


if __name__ == '__main__':
    # parse args
    arg_parser = argparse.ArgumentParser(description="Double strand break repair analyzer (DSBRA)")

    arg_parser.add_argument("-o", "--out_dir", required=True,
                            help="output directory")
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
    arg_parser.add_argument("--full_table", action="store_true",
                            help="Will generate a large table with"
                                 "repair pattern for each analyzed")
#    arg_parser.add_argument( "--intermediate_files", action="store_true",
#                            help="Whether to keep the intermediate files generated by the script")

    args = arg_parser.parse_args()
    print("")

    # set global values
    break_index = args.break_index
    margin = args.margin
    fastq_name = args.fastq
    ref_fa = args.ref_fa
    run_name = fastq_name.split("/")[-1].split(".")[0] \
               + '_b%s_m%s'%(break_index, margin)
    sam_filename = fastq_name[:fastq_name.rfind('.')] + '.sam'
    q_repair_margin = args.q_repair_margin
    last_margin = args.last_margin

    output_dir = args.out_dir
    if os.path.exists(output_dir) :
        #print("Error! The output directory %s already exists. Please use a"
        #      "different name or delete the exisiting output directory."%(output_dir))
        #sys.exit(1)
        print("WARNING! Overwritting the existing output directory %s."%(output_dir))
    else:
        os.makedirs(output_dir)
    print("")

    mut_bamfile_name = os.path.join(output_dir, "%s_mut_reads.sam"%(run_name))
    non_mut_samfile_name = os.path.join(output_dir, "%s_non_mut_reads.sam"%(run_name))
    summary_table_fn = os.path.join(output_dir, "%s_repair_pattern_summary_table.txt"%(run_name))
    if args.full_table:
        repair_pattern_table_fn = os.path.join(output_dir, \
                                               "%s_repair_pattern_full_table.txt"%(run_name))

    output_filename = os.path.join(output_dir, "%s_repair_pattern_summary.txt"%run_name)
    run_info = os.path.join(output_dir, '%s_run_info.txt'%run_name)

    with open(run_info, 'a+') as fh:
        fh.write('Time: %s\nWorking Directory:\n%s\nCommand:\n%s\n'\
                  %(str(datetime.now())[:-7],
                    os.getcwd(),
                    ' '.join(sys.argv)))


    # read in reference fasta, write to new file after checking filetype
    ref_seq = get_ref(ref_fa)

    # store run info
    run_metadata = [('timestamp', str(datetime.now())),
                    ('ref_filename', ref_fa),\
                    ('fastq_file', fastq_name),\
                    ('ref_seq', ref_seq),\
                    ('break_index', break_index), \
                    ('margin', margin)]
    run_metadata = OrderedDict(run_metadata)  # to print in a specific order


    # customed aligner
    run_metadata = aln(sam_filename, ref_fa, fastq_name, break_index, run_metadata,
                       output_dir, record_failed_alignment=True)
    print("")

    # get reads that have mutation events within range
    if os.path.exists(mut_bamfile_name):
        print("Already filtered for reads with mutations within range in %s"%mut_bamfile_name)
    else:
        run_metadata = mut_within_margin(sam_filename, mut_bamfile_name, non_mut_samfile_name,
                          break_index, margin, run_metadata)

    # analyzing the repair pattern of each read, store in a table
    bamfile = pysam.AlignmentFile(mut_bamfile_name+".sorted.bam","rb")

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
        rp_entries = get_rp(read, ref_seq, rp_entries, break_index, margin)

    df = pd.DataFrame(rp_entries, columns=fields)
    if args.full_table:
        print("Writing full table to %s"%repair_pattern_table_fn)
        df.to_csv(repair_pattern_table_fn, sep="\t", index=False)

    # generate summary table of reads. Groups reads with the same repair pattern
    df = get_summary_table(repair_pattern_table_fn,
                           run_metadata, non_mut_samfile_name)
    # the reads made to the end (repair pattern analyzed reads) can be different
    # from reads_with_mutations_within_range
    # because some reads can only mismatches within range, but not immediately close to break index and got filtered out
    run_metadata['Repair pattern analyzed reads'] = df['count'].sum()

    # apply count cutoff
    if args.count_cutoff:
        print("\nApply count cutoff...")
        df = apply_count_cutoff(df, run_metadata,
                                p=args.count_cutoff[0],
                                num_colony=args.count_cutoff[1])

    print("Writing summary table to %s \n"%summary_table_fn)
    df.to_csv(summary_table_fn, sep='\t', index=False)

    # write run metadata
    print("Writing run metadata to %s\n"%(run_info))
    with open(run_info, 'a') as fh:
        for key, item in run_metadata.iteritems():
            if key != 'timestamp':
                fh.write("%s:\n%s\n"%(key, item))

    print("Complete!")
