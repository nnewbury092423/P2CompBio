from math import log2
import csv
from types import AsyncGeneratorType
from threeway_align.BLOSUM62 import *

def compute_gap_penalty(indel_rate,P):
    '''
    This function computes the gap penalty for threeway_align using the BLOSUM62 model and the indel rate
    The transition probabilities have been loaded for you and stored in the nested dictionary P (see file threeway_align/BLOSUM62 and threeway_align/blosum62sym.csv ) 
    '''


    #gap


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

    # arrow is for backtrack
    # initialize (S[i][j][k] = (score,arrow))
    if VERBOSE:
        from sys import stderr; print("Initializing cube", file=stderr)
    S = [[[[0,None] for k in range(len(s3)+1)] for j in range(len(s2)+1)] for i in range(len(s1)+1)]

    S[0][0][0] = [0, None]        # corner
    #print(type(S))
    #print(B)
    # fill in cube axes
    if VERBOSE:
        print("Filling cube axes", file=stderr)
    
    for i in range(1, len(s1)+1): # s1 axis
        #import pdb; pdb.set_trace()
        S[i][0][0][0] = S[i-1][0][0][0] + 2*gap
        S[i][0][0][1] = [1,0,0]
        #Shifted by one
        #print(B)# TODO: replace with your code
    for j in range(1, len(s2)+1): # s2 axis
        S[0][j][0][0] =  S[0][j-1][0][0] + 2*gap 
        S[0][j][0][1] = [0,1,0]
        # TODO: replace with your code
    for k in range(1, len(s3)+1): # s3 axis
        S[0][0][k][0] = S[0][0][k-1][0] + 2*gap # TODO: replace with your code
        S[0][0][k][1] = [0,0,1]
    # fill in cube faces
    if VERBOSE:
        print("Filling cube faces", file=stderr)
    # only 3 faces
    # TODO: add your code here
    # fill 

    for i in range(1, len(s1)+1):
        for j in range(1, len(s2)+1):
            #import pdb; pdb.set_trace()
            option = [S[i-1][j-1][0][0] + B[s1[i-1]][s2[j-1]], S[i][j-1][0][0] + 2*gap, S[i-1][j][0][0] + 2*gap]
            S[i][j][0][0] = max(option)
            max_index = max(enumerate(option), key=lambda x: x[1])[0]
            
            if max_index == 0: 
                arrow = [1, 1, 0]
            elif max_index == 1:
                arrow = [0, 1, 0]
            elif max_index == 2:
                arrow = [1,0, 0]
            S[i][j][0][1] = arrow 



    #import pdb; pdb.set_trace()

    for i in range(1, len(s1)+1):
        for k in range(1, len(s3)+1):
            
            option = [S[i-1][0][k-1][0] + B[s1[i-1]][s3[k-1]], S[i][0][k-1][0] + 2*gap, S[i-1][0][k][0] + 2*gap]
            S[i][0][k][0] = max(option)
            max_index = max(enumerate(option), key=lambda x: x[1])[0]
            
            if max_index == 0: 
                arrow = [1, 0, 1]
            elif max_index == 1:
                arrow = [0,0, 1]
            elif max_index == 2:
                arrow = [1,0,0]

            S[i][0][k][1] = arrow 

    for j in range(1, len(s2)+1):
        for k in range(1, len(s3)+1):

            option = [S[0][j-1][k-1][0] + B[s2[j-1]][s3[k-1]], S[0][j][k-1][0] + 2*gap, S[0][j-1][k][0] + 2*gap]
            S[0][j][k][0] = max(option)
            max_index = max(enumerate(option), key=lambda x: x[1])[0]

            if max_index == 0: 
                arrow = [0, 1, 1]
            elif max_index == 1:
                arrow = [0,0, 1]
            elif max_index == 2:
                arrow = [0,1,0]

            S[0][j][k][1] = arrow 




    # fill in rest of cube
    if VERBOSE:
        print("Filling rest of cube", file=stderr)
    # TODO: add your code here

    for i in range (1, len(s1)+1):   
        for j in range(1, len(s2)+1):
            for k in range(1, len(s3)+1):


                option = [S[i-1][j-1][k-1][0] + B[s1[i-1]][s2[j-1]] + B[s1[i-1]][s3[k-1]] + B[s2[j-1]][s3[k-1]], 
                                 S[i][j-1][k-1][0] + 2*gap + B[s2[j-1]][s3[k-1]],
                                 S[i-1][j][k-1][0] + 2*gap + B[s1[i-1]][s3[k-1]],
                                 S[i-1][j-1][k][0] + 2*gap + B[s1[i-1]][s2[j-1]],
                                 S[i][j][k-1][0] + 2*gap,
                                 S[i][j-1][k][0] + 2*gap,
                                 S[i-1][j][k][0] + 2*gap]

                S[i][j][k][0] = max(S[i-1][j-1][k-1][0] + B[s1[i-1]][s2[j-1]] + B[s1[i-1]][s3[k-1]] + B[s2[j-1]][s3[k-1]], 
                                 S[i][j-1][k-1][0] + 2*gap + B[s2[j-1]][s3[k-1]],
                                 S[i-1][j][k-1][0] + 2*gap + B[s1[i-1]][s3[k-1]],
                                 S[i-1][j-1][k][0] + 2*gap + B[s1[i-1]][s2[j-1]],
                                 S[i][j][k-1][0] + 2*gap,
                                 S[i][j-1][k][0] + 2*gap,
                                 S[i-1][j][k][0] + 2*gap)
                max_index = max(enumerate(option), key=lambda x: x[1])[0]


                if max_index == 0: 
                    arrow = [1, 1, 1]
                elif max_index == 1:
                    arrow = [0,1, 1]
                elif max_index == 2:
                    arrow = [1,0,1]
                elif max_index == 3:
                    arrow = [1,1,0]
                elif max_index == 4:
                    arrow = [0,0,1]
                elif max_index == 5:
                    arrow = [0,1,0]
                elif max_index == 6:
                    arrow = [1,0,0]

                S[i][j][k][1] = arrow 






    
    #bottom back is score
    # backtrack to get alignments
    if VERBOSE:
        print("Backtracking to build alignment", file=stderr)
    
   
        
    #import pdb; pdb.set_trace()
    aln_s1 = ""; aln_s2 = ""; aln_s3 = ""
    arrow = S[-1][-1][-1][1]
    score = S[-1][-1][-1][0]
    #import pdb; pdb.set_trace()
    #print(arrow)
    currentposition = [len(S[0][:][:])-1,len(S[0][:][:])-1,len(S[0][:][:])-1]
    #import pdb; pdb.set_trace()
    while( arrow != None):
        #import pdb; pdb.set_trace()
        if arrow == [0,0,1]:
            aln_s1 = aln_s1 + "_"
            aln_s2 = aln_s2 + "_"
            aln_s3 = aln_s3 + s3[currentposition[2]-1]
        elif arrow == [0,1,0]:
            aln_s1 = aln_s1 + "_"
            aln_s2 = aln_s2 + s2[currentposition[1]-1]
            aln_s3 = aln_s3 + "_"
        elif arrow == [1,0,0]:
            aln_s1 = aln_s1 + s1[currentposition[0]-1]
            aln_s2 = aln_s2 + "_"
            aln_s3 = aln_s3 + "_"
        elif arrow ==  [1,1,0]:
            aln_s1 = aln_s1 + s1[currentposition[0]-1]
            aln_s2 = aln_s2 + s2[currentposition[1]-1]
            aln_s3 = aln_s3 + "_"
        elif arrow == [1,0,1]:
            aln_s1 = aln_s1 + s1[currentposition[0]-1]
            aln_s2 = aln_s2 + "_"
            aln_s3 = aln_s3 + s3[currentposition[2]-1]
        elif arrow ==  [0,1,1]:
            aln_s1 = aln_s1 + "_"
            aln_s2 = aln_s2 + s2[currentposition[1]-1]
            aln_s3 = aln_s3 + s3[currentposition[2]-1]
        elif arrow ==  [1,1,1]:
            #import pdb; pdb.set_trace()
            #print("here")
            aln_s1 = aln_s1 + s1[currentposition[0]-1]
            #import pdb; pdb.set_trace()
            #print("here")
            aln_s2 = aln_s2 + s2[currentposition[1]-1]
            aln_s3 = aln_s3 + s3[currentposition[2]-1]
        
        currentposition = [a - b for a, b in zip(currentposition, arrow)]
        arrow = S[currentposition[0]][currentposition[1]][currentposition[2]][1]
        #print(arrow)

    import pdb; pdb.set_trace()
    #print(score)
    
    return aln_s1[::-1],aln_s2[::-1],aln_s3[::-1],score 

def threeway_align_indel_rate(s1, s2, s3,B,indel_rate,VERBOSE=False):
    gap = compute_gap_penalty(indel_rate,P)                            
    a1,a2,a3,score = threeway_align(s1, s2, s3,B,gap, VERBOSE=VERBOSE) # EXTRA CREDITS: replace with your own algorithm to do alignment
    return a1,a2,a3,score  # optional (extra credits): replace with your own way to do alignment
