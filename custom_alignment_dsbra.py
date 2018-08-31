#!/usr/bin/env python
import re
import sys
import math
import warnings
import argparse
from Bio import pairwise2
from Bio import SeqIO
from Bio.pairwise2 import format_alignment
from itertools import groupby
from joblib import Parallel, delayed
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")

# break site is 0-indexed, the position of the nucleotide right before the break site
# edge distance:
# 1. If the mutation event contains the break site, then edge distance is set to 0
# 2. If the mutation event does not contain the break site, sum of start, end to break site (singed).
# Choose the aln with smallest edge distance. If edge distance is the same, choose the one with smaller start site.
# This essentailly choose the alignment that is cloest to the cut-site as is most biologically plausible.

def aln2cigar(aln):
    ref = aln[0]
    q = aln[1]

    cigar = []
    for i in range(len(ref)):
        if ref[i] == q[i]:
            cigar.append('=')
        else:
            if ref[i] == '-':
                cigar.append('I')
            elif q[i] == '-':
                cigar.append('D')
            else:
                cigar.append('X')

    return cigar


def cal_edge_distance(aln, ref_seq, break_index):
    ref_aln = aln[0]
    query_aln = aln[1]

    if '-' in ref_aln:
        insertions = re.finditer('(-+)', ref_aln)
        for insertion in insertions:
            if insertion.start() <= break_index:
                # update break_index, edge_distance
                break_index += (insertion.end() - insertion.start())

    cigar = aln2cigar(aln)
    muts = list(re.finditer('([^=]+)', ''.join(cigar)))

    for i, m in enumerate(muts):
        if i == 0:  # initate the values
            start = m.start()
            edge_distance = abs(break_index - m.start() + break_index - m.end())

        #if (m.start() <= break_index) and (m.end() >= break_index):
        if  (break_index in range(m.start(), m.end())) or \
            ((break_index + 1) in range(m.start(), m.end())):
            # if the mutation event contains the break site, (break, break+1)
            # then edge distance is set to 0
            tmp_d =0
        else:
            # if the mutation event does not contain the break site,
            # the edge distance is the absolute of the sum of
            # the signed distance of start to break site
            # and the signed distance of end to break site
            tmp_d_1 = abs(break_index - m.start() + break_index - m.end())
            tmp_d_2 = abs(break_index + 1 - m.start() + break_index + 1 - m.end())
            tmp_d = min(tmp_d_1, tmp_d_2)

        # i=0 still goes through this for a correct value of edge_distance
        if tmp_d <= edge_distance:
            edge_distance = tmp_d
            start = m.start()

    return cigar, start, edge_distance


def choose_aln(alns, ref_seq, break_index):
    '''
    choose the aln that is closest to the break site
    '''

    for i, aln in enumerate(alns):
        tmp_cigar, tmp_start, tmp_edge_distance =\
               cal_edge_distance(aln, ref_seq, break_index)

        # initiate the values
        if i == 0:
            cigar = tmp_cigar
            start = tmp_start
            candidate_aln = aln
            edge_distance = tmp_edge_distance

        else:
            if tmp_edge_distance < edge_distance:
                candidate_aln = alns[i]
                start = tmp_start
                cigar = tmp_cigar
                edge_distance  = tmp_edge_distance
            elif tmp_edge_distance == edge_distance:
                if tmp_start < start:
                    candidate_aln = alns[i]
                    start = tmp_start
                    cigar = tmp_cigar

    return candidate_aln, cigar


def trim_aln(aln):
    '''
    For local alignment, ignore anything that extends beyond the ref seq.
    '''
    ref_end = [m.start(0) for m in re.finditer('-*$', aln[0])][0]
    aln = (aln[0][:ref_end], aln[1][:ref_end], aln[2], aln[3], aln[4])
    return aln


def align(X_seq, Y_seq, break_index, idx, align_mode, score_min=(20, 8),
          match=2, mismatch=-6, gap_open=-5, gap_extension=-2):
    '''
    align record X and record Y, generate a SAM format entry
    Args:
        X_seq: reference sequenc, biopython sequence record seq attribute
        Y_seq: query sequence, biopython sequence record seq attribute
        break_index: the position of break
        idx: the record index. It's for keep tracking of how many records the program has processed
        score_min: score function, y = a + b * ln(ref_seq)
        match: match bonus
        mismatch: mismatch penalty
        gap_open: gap open penalty
        gap_extension: gap extension penalty

    '''
    if idx % 500 == 0: print("Aligning seq %s"%idx)

    score_min = score_min[0] + score_min[1] * math.log(len(X_seq))

    if align_mode == "global":
        fw_alns = pairwise2.align.globalms(X_seq, Y_seq, match,
                                           mismatch, gap_open, gap_extension)
        rc_alns = pairwise2.align.globalms(X_seq, Y_seq.reverse_complement(),
                                      match, mismatch, gap_open, gap_extension)
    else:
        fw_alns = pairwise2.align.localms(X_seq, Y_seq, match,
                                           mismatch, gap_open, gap_extension)
        rc_alns = pairwise2.align.localms(X_seq, Y_seq.reverse_complement(),
                                      match, mismatch, gap_open, gap_extension)
        fw_alns = [trim_aln(aln) for aln in fw_alns]
        rc_alns = [trim_aln(aln) for aln in rc_alns]

    # check if the sequence is aligned better with fw or rc
    if fw_alns[0][2] >= rc_alns[0][2]:
        alns = fw_alns
        flag = "0"
        seq = str(Y_seq)
    else:
        alns = rc_alns
        flag = "16"
        seq = str(Y_seq.reverse_complement())

    aln, cigar = choose_aln(alns, X_seq, break_index)
    if align_mode == "local":
        # in local alignemnt the cigar and the query seq may be different length because of the truncation
        seq = ''.join(re.findall('[AGCTN]+', aln[1]))

    #consumed_bases = re.finditer('([^DNHP]+)', ''.join(cigar))
    #pos = str(consumed_bases.next().start() + 1)  # 1-based index
    pos = "1"

    # list cigar to string cigar
    groups = groupby(cigar)
    cigar = ''.join(["%s%s"%(sum(1 for _ in group), label) for label, group in groups])

    score = aln[2]
    if score < score_min:
        flag = "4"
        pos = "0"

    return Y_seq, flag, pos, cigar, seq, score


def get_sam_entry(X, Y, idx, aln_d):
    '''
    align record X and record Y, generate a SAM format entry
    Args:
        X: reference sequence in Biopython sequence record format
        Y: query sequence in Biopython sequence record format
        break_index: the position of break
        idx: the record index. It's for keep tracking of how many records the program has processed
    '''
    if idx%500 == 0:
        print("Record: %s"%idx)

    rname = X.name
    qname = Y.name
    mapq = "255"  # mapq = 255 means mapq is not available
    rnext = "*"
    pnext = "0"
    #tlen = str(len(X.seq))
    tlen = "0"
    qual = "*"

    flag, pos, cigar, seq, score = aln_d[str(Y.seq)]

    sam_entry = "\t".join([qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual]) + '\n'

    return (sam_entry, score)


def main():
    # parse args
    arg_parser = argparse.ArgumentParser(description="Customed aligment for dsbra")
    arg_parser.add_argument("-o", "--outfn", required=True,
                            help="output SAM filename")
    arg_parser.add_argument("--align_mode", required=True, choices=["global",
                            "local"], help="choose the method for alignment")
    arg_parser.add_argument("-f", "--fastq", required=True,
                            help="reads in fastq format for alignment")
    arg_parser.add_argument("--ref", required=True,
                            help="reference sequence")
    arg_parser.add_argument("-b", "--break_site", required=True, type=int,
                            help="The position of break site. 0-based index,\
                                  the postion of nucleotide right before the break site")

    args = arg_parser.parse_args()

    # set break site index
    break_index = args.break_site

    # set output filename
    outfn = args.outfn

    # alignment method
    align_mode = args.align_mode

    # read in reference sequneces and reads
    ref = list(SeqIO.parse(args.ref, "fasta"))[0]
    records = SeqIO.parse(args.fastq, "fastq")
    print("Reference sequence and reads loaded...")

    # set up score_min for alignment
    score_min = (20, 8)

    # uniq_seqs from the records, only perform alignment once for each seq
    # dynamic programming? multiple processing?
    print("Getting the alignment of uniq sequences...")
    uniq_seqs = set()
    for i, r in enumerate(records):
        if i%500000 == 0: print(i)
        uniq_seqs.add(r.seq)
    uniq_seqs = list(uniq_seqs)

    aln_results = Parallel(n_jobs=-1)(delayed(align)\
                                     (ref.seq, uniq_seq, break_index,
                                      idx, align_mode, score_min)
                                      for idx, uniq_seq in enumerate(uniq_seqs))
    aln_d = {}
    for aln_result in aln_results:
        aln_d[str(aln_result[0])] = aln_result[1:]

    # perform the customed alignment and write to a SAM format output file
    with open(outfn, "w+") as outfh:
        print("Writing sam entry for %s reads..."%i)
        # write header
        outfh.write("@HD\tVN:1.0\tSO:unsorted\n")
        outfh.write("@SQ\tSN:%s\tLN:%s\n"%(ref.name, len(ref.seq)))
        outfh.write('''@PG\tID:custom_alignment_dsbra.py\tPN:custom_alignment_dsbra.py\tVN:NA\tCL:"%s"\n'''%(' '.join(sys.argv)))
        # read in the records generator again
        records = SeqIO.parse(args.fastq, "fastq")
        sam_entries_n_score = Parallel(n_jobs=-1, prefer="threads")\
                (delayed(get_sam_entry)(ref, record, idx, aln_d)\
                 for idx, record in enumerate(records))
        sam_entries, scores = map(list, zip(*sam_entries_n_score))
        print("Alignment complete. Writing to %s"%outfn)
        outfh.writelines(sam_entries)
        # plot alignment score distribution
        outplot = outfn + "_alignment_scores.png"
        print("Plotting the alignment score distribution%s"%outplot)
        plt.hist(scores, bins="auto")
        plt.axvline(x = (score_min[0] + score_min[1] * math.log(len(ref.seq))),
                    color="red")
        plt.savefig(outplot)


if __name__ == '__main__':
    main()
