################################################
# 2018-07-31 Ling Chen
#usage: sam_dsbra_interface_ling.py [-h] -o OUT_DIR [--run_name RUN_NAME]
#                                   [-r REF_FA] [-q FASTQ] [-b BREAK_INDEX]
#                                   [-m MARGIN] [--full_table]
#
#Double strand break repair analyzer (DSBRA)
#optional arguments:
#    -h, --help            show this help message and exit
#    -o OUT_DIR, --out_dir OUT_DIR
#                          output directory
#    --run_name RUN_NAME   will be the prefix of output file names
#    -r REF_FA, --ref_fa REF_FA
#                          path of the reference sequence fasta file
#    -q FASTQ, --fastq FASTQ
#                          path of the reads to analyze in fastq format
#    -b BREAK_INDEX, --break_index BREAK_INDEX
#                          break index. 0-indexed, the position of nucleotide
#                          right before the break site.
#    -m MARGIN, --margin MARGIN
#                          the margin before and after the break index that
#                          should be considered for analysis
#    --full_table          Will generate a large table with repair pattern for
#                          each analyzed
################################################
import sys
import argparse
from itertools import groupby
import re
import subprocess
import pysam
import cPickle
import os
import tempfile
import shutil
import numpy as np
import csv
from scipy.stats import norm
import time
from datetime import datetime
from collections import OrderedDict
import pandas as pd
from joblib import Parallel, delayed

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

        self.mut_bottom_range = reference_end  # for the first comparison
        self.mut_top_range = -1  # for the first comparison

        self.mismatch_locations = [ ]
        self.new_mismatch_bases = [ ]
        self.old_mismatch_bases = [ ]
        self.num_transitions = 0
        self.num_transversions = 0

        self.num_insertion_events = 0
        self.insertion_locations = [ ]
        self.inserted_sequences = [ ]

        self.num_deletions = 0
        self.deletion_locations = [ ]
        self.deletion_lengths = [ ]
        self.deleted_sequences = [ ]
        self.deletion_is_microhomologous = [ ]
        self.deletion_micro = [ ]

    def add_insertion(self, ins_loc, ins_str):
        self.num_insertion_events += 1
        self.insertion_locations.append(ins_loc)
        self.inserted_sequences.append(ins_str)

    def add_deletion(self, del_loc, del_size, del_str, is_micro, micro_seq):
        self.num_deletions += 1
        self.deletion_locations.append(del_loc)
        self.deletion_lengths.append(del_size)
        self.deleted_sequences.append(del_str)
        self.deletion_is_microhomologous.append(is_micro)
        self.deletion_micro.append(micro_seq)

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

        self.repair_sequence = repair_seq


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


def mut_within_margin(sam_filename, mut_bamfile_name, non_mut_bamfile_name,
                      break_index, margin, run_metadata):
    '''
    Reads with mutation within margin.
    Args:
        sam_filename: alignment file
        mut_bamfile_name: mutated reads bamfile name
        non_mut_bamfile_name: not mutated reads bamfile name
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
    mut_reads_list = [ ]

    if os.path.exists(mut_bamfile_name):
        print "ERROR: File \"%s\" exists."%mut_bamfile_name
        sys.exit()

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

    run_metadata['coverage'] = int(np.mean(coverage_list))
    run_metadata['num_reads_with_mutations_within_range'] = len(mut_reads_list)

    #Close and sort/index new bamfile
    mut_bamfile.close()
    with open(os.devnull, 'wb') as out:
        subprocess.call(SAMTOOLS_FOLDER+'samtools view -@ 8 -bS '+mut_bamfile_name+' | '+SAMTOOLS_FOLDER+'samtools sort -@ 8 -o '\
            +mut_bamfile_name+'.sorted.bam', shell=True, stdout = out)
        subprocess.call(SAMTOOLS_FOLDER+'samtools index '+ mut_bamfile_name + '.sorted.bam ' + mut_bamfile_name + '.sorted.bai',\
            shell=True, stdout = out)
    bamfile.close()

    run_metadata['num_reads_without_mutations_within_range'] = \
            int(run_metadata['num_success_alignments']) \
            - run_metadata['num_reads_with_mutations_within_range']

    print "First scan and sorting ended. Total time: " + str(time.time()-start_time) + "\n..."

    return run_metadata


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
    is_mmej = False
    microhomology = ''

    matched = True
    while matched:
        before_del_span = (del_start_loc - m, del_start_loc)
        after_del_span = (del_start_loc + del_len, del_start_loc + del_len + m)
        if ref_seq[before_del_span[0] : before_del_span[1]] \
           == ref_seq[after_del_span[0] : after_del_span[1]]:
            m += 1
            is_mmej = True
            mircrohomogy = ref_seq[before_del_span[0] : before_del_span[1]]
        else:
            matched = False

    return is_mmej, microhomology


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
    mut_ranges = []  # based on global locs
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
                    mut_ranges.extend([loc, loc + event_len])
                    for l in range(loc, loc + event_len):
                        rp.add_mismatch(loc, read.query_sequence[q_aln_locs[l]],
                                        ref_seq[r_aln_locs[l]])

            elif event_type == 1:    # insertion
                mut_ranges.extend([loc, loc + event_len])
                ins_str = read.query_sequence[q_loc : q_loc + event_len]
                if loc > 0:  # if the insertion does not occur before ref starts
                    rp.add_insertion(r_aln_locs[loc - 1], ins_str)

            elif event_type == 2:  # deletion
                mut_ranges.extend([loc, loc + event_len])
                del_str = ref_seq[r_loc : r_loc + event_len]
                is_mmej, microhomology = check_microhomology(ref_seq, r_loc, event_len)
                rp.add_deletion(r_loc, event_len, del_str, is_mmej, microhomology)

            else:
                print("Cigar type %s is not implemented. Skip this read: %s"%(event_type, cigar_))
                return "skipped"

        loc += event_len  # update current loc

    mut_ranges = [i for i in mut_ranges if i]
    if mut_ranges:
        rp.mut_bottom_range = min(mut_ranges)
        rp.mut_top_range = max(mut_ranges)
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

            rp_entries['read_name'].append(rp.read_name)
            rp_entries['repair_size'].append(rp.repair_size)
            rp_entries['mismatch_locations'].append(rp.mismatch_locations)
            rp_entries['new_mismatch_bases'].append(rp.new_mismatch_bases)
            rp_entries['num_mismatch'].append((rp.num_transitions + rp.num_transversions))
            rp_entries['num_transitions'].append(rp.num_transitions)
            rp_entries['num_transversions'].append(rp.num_transversions)
            rp_entries['num_insertions'].append(rp.num_insertion_events)
            rp_entries['insertion_locs'].append(rp.insertion_locations)
            rp_entries['insertion_seqs'].append(rp.inserted_sequences)
            rp_entries['num_deletions'].append(rp.num_deletions)
            rp_entries['deletion_locs'].append(rp.deletion_locations)
            rp_entries['deletion_lens'].append(rp.deletion_lengths)
            rp_entries['deletion_sequences'].append(rp.deleted_sequences)
            rp_entries['deletion_is_micro'].append(rp.deletion_is_microhomologous)
            rp_entries['repair_sequence'].append(rp.repair_sequence)
            rp_entries['micro_seq'].append(rp.deletion_micro)

        return rp_entries


if __name__ == '__main__':
    # parse args
    arg_parser = argparse.ArgumentParser(description="Double strand break repair analyzer (DSBRA)")

    arg_parser.add_argument("-o", "--out_dir", required=True,
                            help="output directory")
    arg_parser.add_argument("--run_name", default="test_run",
                            help="will be the prefix of output file names")
    arg_parser.add_argument("-r", "--ref_fa",
                            help="path of the reference sequence fasta file")
    arg_parser.add_argument("-q", "--fastq",
                            help="path of the reads to analyze in fastq format")
    arg_parser.add_argument("-b", "--break_index", type=int,
                            help="break index. 0-indexed, the position of nucleotide right before the break site.")
    arg_parser.add_argument("-m", "--margin", type=int,
                            help="the margin before and after the break index that should be considered for analysis")
    arg_parser.add_argument("--full_table", action="store_true",
                            help="Will generate a large table with repair pattern for each analyzed")
#    arg_parser.add_argument( "--intermediate_files", action="store_true",
#                            help="Whether to keep the intermediate files generated by the script")

    args = arg_parser.parse_args()
    print("")

    # set global values
    break_index = args.break_index
    margin = args.margin
    fastq_name = args.fastq
    ref_fa = args.ref_fa
    run_name = args.run_name
    sam_filename = fastq_name[:fastq_name.rfind('.')]+'.sam'

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
    non_mut_bamfile_name = os.path.join(output_dir, "%s_non_mut_reads.sam"%(run_name))
    summary_table_fn = os.path.join(output_dir, "%s_repair_pattern_summary_table.txt"%(run_name))
    if args.full_table:
        repair_pattern_table_fn = os.path.join(output_dir, \
                                               "%s_repair_pattern_full_table.txt"%(run_name))

    output_filename = os.path.join(output_dir, "%s_repair_pattern_summary.txt"%run_name)
    run_info = os.path.join(output_dir, '%s_run_info.txt'%run_name)

    with open(run_info, 'w+') as fh:
        fh.write('Time: %s\nWorking Directory:\n%s\nCommand:\n%s\n'\
                  %(str(datetime.now())[:-7],
                    os.getcwd(),
                    ' '.join(sys.argv)))


    # read in reference fasta, write to new file after checking filetype
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
        mut_within_margin(sam_filename, mut_bamfile_name, non_mut_bamfile_name,
                          break_index, margin, run_metadata)

    # analyzing the repair pattern of each read, store in a table
    bamfile = pysam.AlignmentFile(mut_bamfile_name+".sorted.bam","rb")

    fields = ['read_name', 'repair_size', 'mismatch_locations',
              'new_mismatch_bases', 'num_mismatch', 'num_transitions',
              'num_transversions', 'num_insertions', 'insertion_locs',
              'insertion_seqs', 'num_deletions', 'deletion_locs',
              'deletion_lens', 'deletion_sequences',
              'deletion_is_micro', 'micro_seq', 'repair_sequence']
    rp_entries = {}
    for f in fields:
        rp_entries[f] = []

    # get_rp is useful is to parallize the code
    for read in bamfile.fetch(bamfile.references[0]):
        rp_entries = get_rp(read, ref_seq, rp_entries, break_index, margin)

    # parallized version
#    rp_entries = Parallel(n_jobs=-1)(delayed(get_rp)(read, ref_seq, rp_entries,
#                                                     break_index, margin)\
#                                     for read in bamfile.fetch(bamfile.references[0]))


    df = pd.DataFrame(rp_entries, columns=fields)
    if args.full_table:
        print("Writing full table to %s"%repair_pattern_table_fn)
        df.to_csv(repair_pattern_table_fn, sep="\t", index=False)

    # generate summary table of reads. Groups reads with the same repair pattern
    print("Writing summary table to %s \n"%summary_table_fn)
    df = pd.read_table(repair_pattern_table_fn)  # this solves the problem caused by a python list in a cell, otherwise it cannot be grouped appropriately
    df = df.drop("read_name", axis=1)
    df['count'] = 1
    df = df.groupby(list(df.columns[:-1]), as_index=False).agg({'count':'sum'})
    df['mut_freq'] = df['count'] / (df['count'].sum())
    df = df.sort_values(by='count', ascending=False)
    # TO DO mut_freq
    df.to_csv(summary_table_fn, sep='\t', index=False)

    # the reads made to the end (repair pattern analyzed reads) can be different
    # from reads_with_mutations_within_range
    # because some reads can only mismatches within range, but not immediately close to break index and got filtered out
    run_metadata['Repair pattern analyzed reads'] = df['count'].sum()

    # write run metadata
    print("Writing run metadata to %s\n"%(run_info))
    with open(run_info, 'a') as fh:
        for key, item in run_metadata.iteritems():
            if key != 'timestamp':
                fh.write("%s:\n%s\n"%(key, item))

    print("Complete!")


###############################################################################
# original aligner uses bowtie2. Fail to get the closet alignment
# to the break site.
###############################################################################


###############################################################################
# original scripts for analyzing the repair pattern of each read
# I guess the logic is:
# for each aligned pairs, if it's a insertion or deletion, use indel to track
# until it finds a matched pair.
# Then it will check under the condition of indel < 0 (deletion),
# whether it's a microhomology or not. (microhomology means nucleotides on
# either end of deletion be the same. mmej)
# It will then check under the condition of indel > 0 (insertion), what the
# insertion str is, and store it. ins_str get reset.
# Finally it will check if the matced pair is a mismatch or not.
# This has potential problem.
# When there is a deletion followed by insertion before a matched pair occur,
# indel might be cancelled out.
###############################################################################
