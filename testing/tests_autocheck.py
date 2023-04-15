import unittest
from threeway_align.todo import *
from threeway_align.utils import *
from os.path import dirname, realpath, join, normpath
from os import listdir, remove
import signal
from subprocess import run
from tempfile import NamedTemporaryFile
import platform

if platform.system() == "Windows":
    from func_timeout import func_timeout, FunctionTimedOut 

TIMEOUT = 180 # seconds

test_path = normpath(join(dirname(realpath(__file__)),"test_data","checking"))
test_cases = ['lowrate','highrate']
test_ids = {'lowrate':['00','01','02','03'],'highrate':['00','01','02','03']}
rate = {'lowrate':0.01, 'highrate':0.1}
cutoff = {'lowrate':0.95,'highrate':0.85}
threshold = 0.8

class Tests(unittest.TestCase):
    def test_01_sanity(self):
        # test reading sample inputs
        for case in test_cases:
            for ID in test_ids[case]:
                test_file = normpath(join(test_path,case,ID + ".fas")) 
            try:
                with open(test_file,'r') as f:
                    seqs = read_FASTA(f)
            except:        
                    self.assertTrue(False,msg="Couldn't read the sample input file " + test_file + "!")

    def test_02_sanity(self):
        # test running FastSP
        for case in test_cases:
            for ID in test_ids[case]:
                test_file = normpath(join(test_path,case,ID + ".aln"))
            tempSP = NamedTemporaryFile(delete=False)
            try:
                run(["java","-jar","FastSP.jar","-r",test_file,"-e",test_file,"-o",tempSP.name],stderr=NamedTemporaryFile(mode='w'),check=True)
                with open(tempSP.name,'r') as f:
                    sp = 0
                    md = 0
                    for line in f:
                        if line.startswith("SP-Score"):
                            sp = float(line.strip().split("SP-Score")[-1])
                        if line.startswith("Modeler"):    
                            md = float(line.strip().split("Modeler")[-1])
                    self.assertTrue(sp == 1.0 and md == 1.0,msg="Wrong run of FastSP on " + test_file)
            except:        
                    self.assertTrue(False,msg="Couldn't run FastSP on " + test_file)

            remove(tempSP.name)

    def __run_case_correctness__(self,ID,gapStr):
        # test if the DP was implemented correctly
        # gap here is either "g1" or "g5" (meaning gap_penalty = -1 or gap_penalty = -5)
        sample_in = normpath(join(test_path, "DP", ID + ".fas"))
        sample_out = normpath(join(test_path, "DP", ID + "." + gapStr + ".aln"))

        self.assertTrue(gapStr in ["g1","g5"],"Invalid gapStr!")
        gap = -1 if gapStr == "g1" else -5
     
        with open(sample_in,'r') as f:
            seqs = read_FASTA(f)
            k1,k2,k3 = sorted(seqs.keys())
        try:
            if platform.system() == "Windows":
                a1,a2,a3,score_student_answer = func_timeout(TIMEOUT,threeway_align,args=(seqs[k1], seqs[k2], seqs[k3],B,gap,))
            else:    
                signal.alarm(TIMEOUT)
                a1,a2,a3,score_student_answer = threeway_align(seqs[k1], seqs[k2], seqs[k3],B,gap)
                signal.alarm(0)
        except: 
            self.assertTrue(False,msg="Failed correctness test on case" + ID + " gap " + gap + ": Couldn't run find_LCAs on input file " + sample_in + "!")
        
        self.assertTrue(len(a1) == len(a2) == len(a3),msg="Failed correctness test on case " + ID + " gap " + str(gap) + ":Aligned sequences does not have equal length")
        for seq,orig in [(a1,seqs[k1]), (a2,seqs[k2]), (a3,seqs[k3])]:
            self.assertTrue(seq.replace('-','') == orig, msg="Failed correctness test on case " + ID + " gap " + str(gap) + ": Removing gaps from the aligned sequences does not yield the original sequences")
        
        with open(sample_out,'r') as f:
            aln = read_FASTA(f)
            a1_true,a2_true,a3_true = aln.values()
            score_true = sum_pairwise_score(a1_true,a2_true,a3_true,gap)    
        self.assertTrue(score_student_answer == score_true,"Failed correctness test on case " + ID + " gap " + str(gap) + ": Wrong pairwise score! Expect: " + str(score_true) + ". Your answer: " + str(score_student_answer))
        score_student_aln = sum_pairwise_score(a1,a2,a3,gap)
        self.assertTrue(score_student_aln == score_true,"Failed correctness test on case " + ID + " gap " + str(gap) + ": Alignment got wrong pairwise score! Expect: " + str(score_true) + ". Your alignment score: " + str(score_student_aln))
        
    def __run_case_accuracy__(self,case,ID):
        sample_in = normpath(join(test_path,case,ID + ".fas"))
        sample_ref = normpath(join(test_path,case, ID + ".aln"))
        
        with open(sample_in,'r') as f:
            seqs = read_FASTA(f)
            k1,k2,k3 = sorted(seqs.keys())
        try:
            if platform.system() == "Windows":
                a1,a2,a3,_ = func_timeout(TIMEOUT,threeway_align_indel_rate,args=(seqs[k1], seqs[k2], seqs[k3], B,rate[case],))
            else:    
                signal.alarm(TIMEOUT)
                a1,a2,a3,_ = threeway_align_indel_rate(seqs[k1], seqs[k2], seqs[k3], B,rate[case])
                signal.alarm(0)
        except: 
            self.assertTrue(False,msg="Failed accuracy test on " + case + " " + ID + ": Couldn't run threeway_align  on input file " + sample_in + "!")
        
        tempOut = NamedTemporaryFile(delete=False)
        
        with open(tempOut.name,'w') as f:
            for k,a in [(k1,a1),(k2,a2),(k3,a3)]:
                f.write(">"+k+"\n"+a+"\n")
        
        tempSP = NamedTemporaryFile(delete=False)
        run(["java","-jar","FastSP.jar","-r",sample_ref,"-e",tempOut.name,"-o",tempSP.name],stderr=NamedTemporaryFile(mode='w'),check=True)
        
        with open(tempSP.name,'r') as fin:
            for line in fin:
                if line.startswith("SP-Score"):
                    score = float(line.split("SP-Score")[-1].strip())
                    break

        self.assertTrue(score >= threshold, msg = "Failed accuracy test on" + case + " " + ID + ": Too low SP-Score. Expect: >=" + str(threshold) + ". Your score: " + str(score))
        if score < cutoff[case]:
            print("Partially passed accuracy test on" + case + " " + ID + ": Your SP-Score is " + str(score) + ". Aim at score " + str(cutoff[case]) + " to get full credit") 
        remove(tempOut.name)
        remove(tempSP.name)
    
    def test_03_correctness(self):
        ID = '00'
        gapStr = 'g1'
        self.__run_case_correctness__(ID,gapStr)

    def test_04_correctness(self):
        ID = '00'
        gapStr = 'g5'
        self.__run_case_correctness__(ID,gapStr)

    def test_05_correctness(self):
        ID = '01'
        gapStr = 'g1'
        self.__run_case_correctness__(ID,gapStr)
    
    def test_06_correctness(self):
        ID = '01'
        gapStr = 'g5'
        self.__run_case_correctness__(ID,gapStr)
    
    def test_07_correctness(self):
        ID = '02'
        gapStr = 'g1'
        self.__run_case_correctness__(ID,gapStr)
    
    def test_08_correctness(self):
        ID = '02'
        gapStr = 'g5'
        self.__run_case_correctness__(ID,gapStr)
    
    def test_09_accuracy(self):
        case = 'lowrate'
        ID = '01'
        self.__run_case_accuracy__(case,ID)
    
    def test_10_accuracy(self):
        case = 'lowrate'
        ID = '02'
        self.__run_case_accuracy__(case,ID)

    def test_11_accuracy(self):
        case = 'lowrate'
        ID = '03'
        self.__run_case_accuracy__(case,ID)
    
    def test_12_accuracy(self):
        case = 'highrate'
        ID = '01'
        self.__run_case_accuracy__(case,ID)
    
    def test_13_accuracy(self):
        case = 'highrate'
        ID = '02'
        self.__run_case_accuracy__(case,ID)

    def test_14_accuracy(self):
        case = 'highrate'
        ID = '03'
        self.__run_case_accuracy__(case,ID)
