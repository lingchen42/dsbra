#!/usr/bin/env python
import re
import sys
import math
import gzip
import warnings
import argparse
from scipy import stats
import pandas as pd
from Bio import pairwise2
from Bio import SeqIO
from Bio.Seq import Seq
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


def repetitions(s):
       r = re.compile(r"(.+?)\1+$")
       for match in r.finditer(s):
            return len(match.group(0))/len(match.group(1))


def remove_pcr_tail(q_seq, pcr_tail_seq, min_len, pcr_primer1, pcr_primer2,
                    check_n_nuc_at_end, min_score=20,
                    max_end_repeats=5, fuzzy_len=8):
    q_seq = str(q_seq)
    # if find a perfect match for the pcr_tail
    m = [m.start() for m in re.finditer(pcr_tail_seq, q_seq)]
    if m:
        # truncate q_seq to only keep anything before the pcr_tail_seq
        q_seq = q_seq[:m[0]]

    # sometimes there is not a perfect match but see a lot of repeats at the end
    # then look more closely for possible pcr_tail_seq
    #elif repetitions(q_seq) > max_end_repeats:

    # sometimes there is not a perfect match but a fuzzy match
    else:
        aln = pairwise2.align.localms(pcr_tail_seq, q_seq, 2, -6, -5, -2)[0]
        if aln[2] > min_score:
            # find the match > fuzzy_len bps
            fuzzym = [m.start() for m in \
             re.finditer("(=){%s,}"%fuzzy_len, ''.join(aln2cigar(aln)))]
            if fuzzym:
                # truncate q_seq to only keep anything before the pcr_tail_seq
                try:
                    q_seq = re.findall("[AGCTN]+", aln[1][:fuzzym[0]])[0]
                except:
                    # in case no match
                    q_seq = ''

    out_seq = False
    if len(q_seq) > min_len:
        # require the read to start with pcr primer seq (fuzzy matched)
        # and end with part of the reverse primer seq
        primer_match_score = pairwise2.align.localms(pcr_primer1,
                                                     q_seq[:len(pcr_primer1)],
                                                     2, -6, -5, -2)[0][2]
        if primer_match_score > min_score:
            if pcr_primer2 != '-1':  # if we are going to check reverse primer
                if (q_seq[-check_n_nuc_at_end:] in pcr_primer2):
                    out_seq = Seq(q_seq)
            else:
                out_seq = Seq(q_seq)

    return out_seq


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


def is_large_del(cigar, max_mismatches=3):
    '''
    large deletion should have cigar string like:
        37=116D18=1I
        or
        25=52D4=20D70=1I
        or
        25=52D4=20D15=1X14=1X39=1I
    criteria 1: The sequence in seq should always be found in ref_seq allowing
              some mismatches, that is n_mismatches < some cutoff.
    criterial 2: Should have part of primer 1 at the begining (already
              filtered during valid/not valid read) part.
    '''
    n_mismatches = sum([int(i) for i in re.findall('([0-9]+)X', cigar)])
    if n_mismatches <= max_mismatches:
        return True
    else:
        return False


def align(X_seq, Y_seq, break_index, idx, align_mode, pcr_tail_seq, min_len,
          pcr_primer1, pcr_primer2, check_n_nuc_at_end,
          score_min=(20, 8), strand="Unknown", match=2, mismatch=-6,
          gap_open=-5, gap_extension=-2):
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
    if idx % 1000 == 0: print("Aligning seq %s"%idx)

    ori_y_seq = Y_seq
    Y_seq = remove_pcr_tail(Y_seq, pcr_tail_seq, min_len,
                            pcr_primer1, pcr_primer2, check_n_nuc_at_end)
    if Y_seq != False:  # if it's a valid read, align it
        score_min = score_min[0] + score_min[1] * math.log(len(X_seq))

        # check if the sequence is aligned better with fw or rc
        if strand == "Unknown":
            if align_mode == "global":
                fw_alns = pairwise2.align.globalms(X_seq, Y_seq, match,
                                               mismatch, gap_open, gap_extension)
                rc_alns = pairwise2.align.globalms(X_seq, Y_seq.reverse_complement(),
                                          match, mismatch, gap_open, gap_extension)
            else:
                fw_alns = pairwise2.align.localms(X_seq, Y_seq, match,
                                               mismatch, gap_open, gap_extension)
                fw_alns = [trim_aln(aln) for aln in fw_alns]
                rc_alns = pairwise2.align.localms(X_seq, Y_seq.reverse_complement(),
                                          match, mismatch, gap_open, gap_extension)
                rc_alns = [trim_aln(aln) for aln in rc_alns]

            if fw_alns[0][2] >= rc_alns[0][2]:
                alns = fw_alns
                flag = "0"
                seq = str(Y_seq)
            else:
                alns = rc_alns
                flag = "16"
                seq = str(Y_seq.reverse_complement())

        elif strand == "forward":
            flag = "0"
            seq = str(Y_seq)
            if align_mode == "global":
                alns = pairwise2.align.globalms(X_seq, Y_seq, match,
                                               mismatch, gap_open, gap_extension)
            else:  # local
                alns = pairwise2.align.localms(X_seq, Y_seq, match,
                                               mismatch, gap_open, gap_extension)
                alns = [trim_aln(aln) for aln in alns]

        else:  # reverse
            flag = "16"
            seq = str(Y_seq.reverse_complement())
            if align_mode == "global":
                alns = pairwise2.align.globalms(X_seq, Y_seq.reverse_complement(),
                                          match, mismatch, gap_open, gap_extension)
            else:  # local
                alns = pairwise2.align.localms(X_seq, Y_seq.reverse_complement(),
                                          match, mismatch, gap_open, gap_extension)
                alns = [trim_aln(aln) for aln in alns]

        aln, cigar = choose_aln(alns, X_seq, break_index)
        if align_mode == "local":
            # in local alignemnt the cigar and the query seq may be different length because of the truncation
            seq = ''.join(re.findall('[AGCTN]+', aln[1]))

        #consumed_bases = re.finditer('([^DNHP]+)', ''.join(cigar))
        #pos = str(consumed_bases.next().start() + 1)  # 1-based index
        pos = "1"

        # list cigar to string cigar
        groups = groupby(cigar)
        cigar_str = ''.join(["%s%s"%(sum(1 for _ in group), label) for label, group in groups])

        score = aln[2]
        if score < score_min:
            if not is_large_del(cigar_str):
                flag = "4"
                pos = "0"

        return ori_y_seq, flag, pos, cigar_str, seq, score

    else:
        return ori_y_seq, None


def align_wrapper(X_seq, Y_seq, break_index, idx, align_mode, pcr_tail_seq, min_len,
          pcr_primer1, pcr_primer2, check_n_nuc_at_end,
          score_min=(20, 8), strand="Unknown", cigar_len_check=30,
          match=2, mismatch=-6, gap_open=-5, gap_extension=-2):
    '''
    Perform check on the cigar string. If too complicated, try use the other
    alignment method
    '''
    aln_results = align(X_seq, Y_seq, break_index,
                        idx, align_mode, pcr_tail_seq, min_len,
                        pcr_primer1, pcr_primer2,
                        check_n_nuc_at_end,
                        score_min, strand)

    opposite_align = {'global':'local', 'local':'global'}
    if aln_results[1] != None:
        cigar = aln_results[3]
        # check if cigar is too complicated
        # if it is, use another alignment method
        if len(cigar) > cigar_len_check:
            align_mode = opposite_align[align_mode]
#            print("complicated cigar detected %s, try %s alignment"%(cigar,
#                                                                     align_mode))
            other_aln_results = align(X_seq, Y_seq, break_index,
                                idx, align_mode, pcr_tail_seq, min_len,
                                pcr_primer1, pcr_primer2,
                                check_n_nuc_at_end,
                                score_min, strand)
            other_cigar = other_aln_results[3]
            if (len(other_cigar) < len(cigar)) and\
               (re.findall('^[0-9]+(.{1})', other_cigar)[0] != 'D'):   # don't start with deletion
#                print('simpler cigar %s found with %s'%(other_cigar, align_mode))
                aln_results = other_aln_results

    return aln_results


def get_sam_entry(X, Y, pcr_tail_seq, min_len, idx, aln_d):
    '''
    align record X and record Y, generate a SAM format entry
    Args:
        X: reference sequence in Biopython sequence record format
        Y: query sequence in Biopython sequence record format
        break_index: the position of break
        idx: the record index. It's for keep tracking of how many records the program has processed
    '''
    if idx%5000 == 0:
        print("Record: %s"%idx)

    rname = X.name
    qname = Y.name
    mapq = "255"  # mapq = 255 means mapq is not available
    rnext = "*"
    pnext = "0"
    #tlen = str(len(X.seq))
    tlen = "0"
    qual = "*"

    aln_result = aln_d[str(Y.seq)]
    if aln_result[0] != None:
        flag, pos, cigar, seq, score = aln_result
        sam_entry = "\t".join([qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual]) + '\n'
        return (sam_entry, score)
    else:
        return None


def main():
    # parse args
    arg_parser = argparse.ArgumentParser(description="Customed aligment for dsbra")
    arg_parser.add_argument("-o", "--outfn", required=True,
                            help="output SAM filename")
    arg_parser.add_argument("--not_valid_reads", default=False,
                            help="output not valid reads")
    arg_parser.add_argument("--run_info", required=True,
                            help="store alignment statistics")
    arg_parser.add_argument("--align_mode", required=True, choices=["global",
                            "local"], help="choose the method for alignment")
    arg_parser.add_argument("-f", "--fastq", required=True,
                            help="reads in fastq format for alignment")
    arg_parser.add_argument("--ref", required=True,
                            help="reference sequence")
    arg_parser.add_argument("--pcr_primer1", required=True,
                            help="primer seq, at the begining of the ref seq; will require the read to starts with this primer seq")
    arg_parser.add_argument("--pcr_primer2", default=-1,
                            help="primer seq, at the end of the ref seq; will require last 4 nucleotide of the read is part of this primer seq; if not used, use -1")
    arg_parser.add_argument("--pcr_tail_seq",
                            default="ATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
                            help="added sequence for PCR;default ATCGGAAGAGCACACGTCTGAACTCCAGTCAC")
    arg_parser.add_argument("--min_len", type=int, default=44,
                            help="miminal valid length (reads - pcr_tail_seq); default 44bp")
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

    # run info
    run_info = args.run_info

    # not valid reads
    not_valid_reads = args.not_valid_reads

    # read in reference sequneces and reads
    ref = list(SeqIO.parse(args.ref, "fasta"))[0]
    if args.fastq[-2:] == 'gz':
        fastq_fh = gzip.open(args.fastq, "rt")
        records = SeqIO.parse(fastq_fh, "fastq")
    else:
        records = SeqIO.parse(args.fastq, "fastq")
    print("Reference sequence and reads loaded...")

    # params for removing pcr tail seq
    pcr_tail_seq = args.pcr_tail_seq
    min_len = args.min_len
    print("Will use %s and minimum length of %sbp to filter out primer/dimer"%(pcr_tail_seq, min_len))

    # params for checking primer existance
    pcr_primer1 = args.pcr_primer1
    pcr_primer2 = args.pcr_primer2
    check_n_nuc_at_end = 4
    print("Will require a fuzzy match of %s at the begining of the seq"%pcr_primer1)
    print("Will require last 4 nucleotides in the read in %s"%pcr_primer2)

    # set up score_min for alignment
    score_min = (20, 8)

    # uniq_seqs from the records, only perform alignment once for each seq
    # dynamic programming? multiple processing?
    print("Getting uniq seqs....")
    uniq_seqs = set()
    for i, r in enumerate(records):
        if i%50000 == 0: print(i)
        uniq_seqs.add(r.seq)
    uniq_seqs = list(uniq_seqs)

    print("\nGetting the alignment of %s uniq reads..."%(len(uniq_seqs)))
    # check the first 100 for strand
    strand = "Unknown"
    seq_to_check_start = 0
    seq_to_check_end = 100
    while strand == "Unknown":  # in case the first 100 doesnot align
        aln_results = Parallel(n_jobs=-1)(delayed(align_wrapper)\
                                         (ref.seq, uniq_seq, break_index,
                                          idx, align_mode, pcr_tail_seq, min_len,
                                          pcr_primer1, pcr_primer2,
                                          check_n_nuc_at_end,
                                          score_min, strand)\
                                          for idx, uniq_seq in\
                     enumerate(uniq_seqs[seq_to_check_start:seq_to_check_end]))
        aln_d = {}
        flags = []
        for aln_result in aln_results:
            aln_d[str(aln_result[0])] = aln_result[1:]
            if aln_result[1] not in [None, "4"]:
                flags.append(aln_result[1])

        if flags:
            if stats.mode(flags) == "16":
                strand = "reverse"
            else:
                strand = "forward"
        else:
            seq_to_check_start += 100
            seq_to_check_end += 100
    print("\nStrand %s detacted using first %s sequences\n"%(strand, seq_to_check_end))

    # do the alignment for the rest of seqs
    aln_results = Parallel(n_jobs=-1)(delayed(align_wrapper)\
                                     (ref.seq, uniq_seq, break_index,
                                      idx, align_mode, pcr_tail_seq, min_len,
                                      pcr_primer1, pcr_primer2,
                                      check_n_nuc_at_end,
                                      score_min, strand)\
                                      for idx, uniq_seq in\
                                      enumerate(uniq_seqs[seq_to_check_end:]))
    if not_valid_reads:
        print("Writing not valid reads to %s"%not_valid_reads)
        with open(not_valid_reads, 'w+') as not_valid_reads_fh:
            for aln_result in aln_results:
                aln_d[str(aln_result[0])] = aln_result[1:]
                if not aln_result[1]:  # output not valid reads
                    not_valid_reads_fh.write(str(aln_result[0]) + '\n')
    else:
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
        sam_entries_n_score = Parallel(n_jobs=-1, backend="threading")\
                (delayed(get_sam_entry)(ref, record, pcr_tail_seq, min_len,
                                        idx, aln_d)\
                 for idx, record in enumerate(records))
        print("Alignment complete. Writing to %s"%outfn)
        sam_entries_n_score = [s for s in sam_entries_n_score if s!=None]
        sam_entries, scores = map(list, zip(*sam_entries_n_score))
        outfh.writelines(sam_entries)

        # save run info
        df_run_info = pd.read_csv(run_info)
        df_run_info["n_primerdimer"] = i - len(sam_entries_n_score)
        df_run_info["n_non_primerdimer_reads"] = len(sam_entries_n_score)
        df_run_info.to_csv(run_info, index=False)

        # plot alignment score distribution
        outplot = outfn + "_alignment_scores.png"
        print("Plotting the alignment score distribution %s"%outplot)
        plt.hist(scores, bins="auto")
        plt.axvline(x = (score_min[0] + score_min[1] * math.log(len(ref.seq))),
                    color="red")
        plt.savefig(outplot)


if __name__ == '__main__':
    main()
