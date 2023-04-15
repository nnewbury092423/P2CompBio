from math import log2
import csv
from threeway_align.BLOSUM62 import *

def compute_gap_penalty(indel_rate,P):
    '''
    This function computes the gap penalty for threeway_align using the BLOSUM62 model and the indel rate
    The transition probabilities have been loaded for you and stored in the nested dictionary P (see file threeway_align/BLOSUM62 and threeway_align/blosum62sym.csv ) 
    '''
    gap = 0 # TODO: replace with your way to compute gap penalty
    return gap

def threeway_align(s1, s2, s3, B, gap, VERBOSE=False):
    '''
    This function computes the threeway sequence alignment between strings s1, s2, and s3 given a substitution matrix and gap penalty
    :param s1 (string): the first sequence
    :param s2 (string): the second sequence
    :param s3 (string): the third sequence
    :param B (char,char -> int) : the substitution matrix (e.g. BLOSUM62)
    :param gap (negative float):  gap penalty
    '''
    # initialize (S[i][j][k] = (score,arrow))
    if VERBOSE:
        from sys import stderr; print("Initializing cube", file=stderr)
    S = [[[None for k in range(len(s3)+1)] for j in range(len(s2)+1)] for i in range(len(s1)+1)]
    S[0][0][0] = (0, None)        # corner

    # fill in cube axes
    if VERBOSE:
        print("Filling cube axes", file=stderr)
    for i in range(1, len(s1)+1): # s1 axis
        S[i][0][0] = 0 # TODO: replace with your code
    for j in range(1, len(s2)+1): # s2 axis
        S[0][j][0] = 0 # TODO: replace with your code
    for k in range(1, len(s3)+1): # s3 axis
        S[0][0][k] = 0 # TODO: replace with your code

    # fill in cube faces
    if VERBOSE:
        print("Filling cube faces", file=stderr)
    # TODO: add your code here

    # fill in rest of cube
    if VERBOSE:
        print("Filling rest of cube", file=stderr)
    # TODO: add your code here

    # backtrack to get alignments
    if VERBOSE:
        print("Backtracking to build alignment", file=stderr)
    aln_s1 = ""; aln_s2 = ""; aln_s3 = ""
 
    score = S[-1][-1][-1]
    return aln_s1[::-1],aln_s2[::-1],aln_s3[::-1],score 

def threeway_align_indel_rate(s1, s2, s3,B,indel_rate,VERBOSE=False):
    gap = compute_gap_penalty(indel_rate,P)                            
    a1,a2,a3,score = threeway_align(s1, s2, s3,B,gap, VERBOSE=VERBOSE) # EXTRA CREDITS: replace with your own algorithm to do alignment
    return a1,a2,a3,score  # optional (extra credits): replace with your own way to do alignment
