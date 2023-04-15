#!/usr/bin/env python3

from threeway_align.todo import *
from threeway_align.utils import *

if __name__ == "__main__":
    '''
    This is how we handle loading the input dataset, running your function, and printing the output
    '''
    import argparse
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input_alignment', required=False, type=str, default='stdin', help="Input Alignment (FASTA format)")
    parser.add_argument('-g', '--gap', required=True, type=float, help="Gap penalty")
    args = parser.parse_args()
    if args.input_alignment == 'stdin':
        from sys import stdin as infile
    else:
        infile = open(args.input_alignment)
    seqs = read_FASTA(infile); infile.close()
    assert len(seqs) == 3, "Input FASTA file must have exactly 3 sequences"
    k1,k2,k3 = sorted(seqs.keys())
    print(sum_pairwise_score(seqs[k1],seqs[k2],seqs[k3],args.gap))
