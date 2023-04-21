#!/usr/bin/env python3

from threeway_align.todo import *
from threeway_align.utils import *
from threeway_align.BLOSUM62 import *

if __name__ == "__main__":
    '''
    This is how we handle loading the input dataset, running your function, and printing the output
    '''
    import argparse
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input_sequences', required=False, type=str, default='stdin', help="Input Sequences (FASTA format)")
    parser.add_argument('-r', '--indel_rate', required=False, type=float, default=0.1, help="Indel Rate")
    parser.add_argument('-g', '--gap_penalty', required=False, type=float, help="Gap Penalty. Must be a negative number. Will override indel rate")
    parser.add_argument('-o', '--output_alignment', required=False, type=str, default='stdout', help="Output Alignment (FASTA format)")
    args = parser.parse_args()
    assert args.indel_rate >= 0, "Relative indel rate must be non-negative"
    if args.gap_penalty is not None:
        assert args.gap_penalty <0 , "Gap penalty must be negative"
    if args.input_sequences == 'stdin':
        from sys import stdin as infile
    else:
        infile = open(args.input_sequences)
    if args.output_alignment == 'stdout':
        from sys import stdout as outfile
    else:
        outfile = open(args.output_alignment,'w')
    seqs = read_FASTA(infile); infile.close()
    assert len(seqs) == 3, "Input FASTA file must have exactly 3 sequences"
    k1,k2,k3 = sorted(seqs.keys())
    if args.gap_penalty is not None:
        a1,a2,a3,score = threeway_align(seqs[k1], seqs[k2], seqs[k3],B,args.gap_penalty)
    else:
        a1,a2,a3,score = threeway_align_indel_rate(seqs[k1], seqs[k2], seqs[k3], B, args.indel_rate)
    assert len(a1) == len(a2) == len(a3), "Aligned sequences must have equal length"
    assert len(a1) != 0, "Returned empty strings"
    for seq,orig in [(a1,seqs[k1]), (a2,seqs[k2]), (a3,seqs[k3])]:
        assert seq.replace('-','') == orig, "Removing gaps from the aligned sequences must yield the original sequences"
    for ID,seq in [(k1,a1), (k2,a2), (k3,a3)]:
        outfile.write(">%s\n%s\n" % (ID,seq))
    print("Optimal pairwise score: " + str(score))
    outfile.close()
