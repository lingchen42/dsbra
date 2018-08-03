#!/usr/bin/env python
import sys
import math
from Bio import pairwise2
from Bio import SeqIO
from Bio.pairwise2 import format_alignment
import re
import argparse
from itertools import groupby
from joblib import Parallel, delayed

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
            if insertion.start() < break_index:
                # update break_index, edge_distance
                break_index += (insertion.end() - insertion.start())

    cigar = aln2cigar(aln)
    muts = re.finditer('([^=]+)', ''.join(cigar))
    edge_distance = len(ref_seq)
    for m in muts:
        if (m.start() < break_index) and (m.end() > break_index):
            # if the mutation event contains the break site, then edge distance is set to 0
            tmp_d =0
        else:
            # if the mutation event does not contain the break site,
            # the edge distance is the absolute of the sum of the signed distance of start to break site
            # and the signed distance of end to break site
            tmp_d = abs(break_index - m.start() + break_index - m.end())

        if tmp_d < edge_distance:
            edge_distance = tmp_d
            start = m.start()

    return cigar, start, edge_distance


def choose_aln(alns, ref_seq, break_index):
    # choose the aln that is closest to the break site
    start = 0
    cigar = []
    edge_distance = len(ref_seq)

    for i, aln in enumerate(alns):
        tmp_cigar, tmp_start, tmp_edge_distance = cal_edge_distance(aln, ref_seq, break_index)
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


def align(X, Y, break_index, idx, score_min=(20, 8), match=2, mismatch=-6, gap_open=-5, gap_extension=-2):
    '''
    align record X and record Y, generate a SAM format entry
    Args:
        X: reference sequence in Biopython sequence record format
        Y: query sequence in Biopython sequence record format
        break_index: the position of break
        idx: the record index. It's for keep tracking of how many records the program has processed
        score_min: score function, y = a + b * ln(ref_seq)
        match: match bonus
        mismatch: mismatch penalty
        gap_open: gap open penalty
        gap_extension: gap extension penalty

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

    score_min = score_min[0] + score_min[1] * math.log(len(X))

    fw_alns = pairwise2.align.globalms(X.seq, Y.seq, match, mismatch, gap_open, gap_extension)
    rc_alns = pairwise2.align.globalms(X.seq, Y.seq.reverse_complement(),
                                  match, mismatch, gap_open, gap_extension)

    # check if the sequence is aligned better with fw or rc
    if fw_alns[0][2] >= rc_alns[0][2]:
        alns = fw_alns
        flag = "0"
        seq = str(Y.seq)
    else:
        alns = rc_alns
        flag = "16"
        seq = str(Y.seq.reverse_complement())

    aln, cigar = choose_aln(alns, X.seq, break_index)
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

    sam_entry = "\t".join([qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual]) + '\n'

    return sam_entry


def main():
    # parse args
    arg_parser = argparse.ArgumentParser(description="Customed aligment for dsbra")
    arg_parser.add_argument("-o", "--outfn", required=True,
                            help="output SAM filename")
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

    # read in reference sequneces and reads
    ref = list(SeqIO.parse(args.ref, "fasta"))[0]
    records = list(SeqIO.parse(args.fastq, "fastq"))
    print("Reference sequence and %s reads loaded. Begin alignment..."%(len(records)))

    # perform the customed alignment and write to a SAM format output file
    with open(outfn, "w+") as outfh:
        # write header
        outfh.write("@HD\tVN:1.0\tSO:unsorted\n")
        outfh.write("@SQ\tSN:%s\tLN:%s\n"%(ref.name, len(ref.seq)))
        outfh.write('''@PG\tID:custom_alignment_dsbra.py\tPN:custom_alignment_dsbra.py\tVN:NA\tCL:"%s"\n'''%(' '.join(sys.argv)))
        sam_entries = Parallel(n_jobs=-1)(delayed(align)(ref, record, break_index, idx) for idx, record in enumerate(records))
        print("Alignment complete. Writing to %s"%outfn)
        outfh.writelines(sam_entries)


if __name__ == '__main__':
    main()
