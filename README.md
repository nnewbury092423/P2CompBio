[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-24ddc0f5d75046c5622901739e7c5dd533143b0c8e959d652212380cedb1ea36.svg)](https://classroom.github.com/a/3LGf16Zh)
# ECE 208 Homework 2: Three-Way Sequence Alignment and Inferring Gap Penalty

## Part 1: Three-way sequence alignment
We extend the pairwise sequence alignment algorithm to three sequences. You will be given three amino acid sequences *s*<sub>1</sub>, *s*<sub>2</sub>, and *s*<sub>3</sub> with sequences lengths *m*<sub>1</sub>, *m*<sub>2</sub>, and *m*<sub>3</sub>, respectively. The sequences are homologous, but have substitutions, insertions, and deletions, and therefore, differ in lengths.

### Problem
Design a **dynamic programming algorithm** to align 3 sequences *s*<sub>1</sub>, *s*<sub>2</sub>, and *s*<sub>3</sub> such that **the sum of pairwise alignment scores** is maximized, given the **BLOSUM62**<sup>1</sup> substitution model and **non-affine gap penalty**.

### Input
The `threeway_align` function takes as input the following parameters:
* Three Python [strings](https://docs.python.org/3/library/stdtypes.html#textseq) `s1`, `s2`, and `s3`, representing the three amino acid sequences *s*<sub>1</sub>, *s*<sub>2</sub>, and *s*<sub>3</sub>
* A real-valued number (usually is negative) `gap`, which is the gap penalty.

### Output
The `threeway_align` function returns three strings `aligned_s1`, `aligned_s2`, and `aligned_s3` with the following properties:
* All three have equal length
    * `len(aligned_s1) == len(aligned_s2) == len(aligned_s3)`
* Removing all gap characters from `aligned_s1`, `aligned_s2`, and `aligned_s3` yields `s1`, `s2`, and `s3`, respectively
    * `aligned_s1.replace('-','') == s1 and aligned_s2.replace('-','') == s2 and aligned_s3.replace('-','') == s3`

### Hints
If you compare this problem to the pairwise alignment problem, you will see that instead of a 2D matrix as in the pairwise case, you will need to compute a 3D matrix of dimensions (|*s*<sub>1</sub>|+1) × (|*s*<sub>2</sub>|+1) × (|*s*<sub>3</sub>|+1). The recursive formula will also be similar, just with more cases involved.

## Part 2: Inferring gap penalty from the substitution model and indel rate
Recall that the indel rate is the probability of an indel relative to the probability of a substitution. For example, an indel rate of 0.1 means that, on average, for every 10 substitutions, you would expect 1 insertion or deletion. 
In this part of the assignment you are given the **indel rate** instead of gap penalty,and you need to find the best settings of gap penalty with respect to the given indel rate and the **BLOSUM62** model.

### Hints
When you want to find gap penalty, you may find that you want to know the *p*<sub>*ij*</sub> values used in the BLOSUM62 matrix, which we have included in [blosum62sym.csv](blosum62sym.csv). The other thing you may need to know is how this is turned into the BLOSUM62 matrix. Convince yourself that what they used is the following (important information here is that log is in base-2 and there is a multiplication by 2):

![equation](https://latex.codecogs.com/svg.latex?2%5Clog_%7B2%7D%5Cleft%28%5Cfrac%7Bp_%7Bij%7D%7D%7Bq_%7Bi%7Dq_%7Bj%7D%7D%5Cright%20%29)

## Starting Code
We will use the [FASTA format](https://en.wikipedia.org/wiki/FASTA_format) for both inputs and outputs. For your convenience, we provide the code for reading and writing FASTA files. Your task is to complete the file [threeway_align/todo.py](threeway_align/todo.py).
The package given to you includes the following components:

* **[threeway_align/todo.py](threeway_align/todo.py):** You will need to address all the places marked with TODO. 
In this assignment, there are two functions you need to complete: `threeway_align` (Part 1) and `compute_gap_penalty` (Part 2)
* **[threeway_align.py](threeway_align.py):** Given to you. Please **DO NOT** modify. This file handles I/O for your program. 
Type `python3 threeway_align.py -h` to see the input/output options.
* **[compute_pairwise_score.py](compute_pairwise_score.py):** Given to you. Please **DO NOT** modify. This script is provided for testing and debugging purposes. The script computes the sum-of-pairwise score of an alignment. Type `python3 compute_pairwise_score.py -h` to see the input/output options.
* **[threeway_align/utils.py](threeway_align/utils.py):** Given to you. Please **DO NOT** modify. 
This file contains the functions `read_FASTA` and `sum_pairwise_score` that will be used to run and test your program. 
* **[threeway_align/BLOSUM62.py](threeway_align/BLOSUM62.py):** Given to you. Please **DO NOT** modify. This file provides the BLOSUM62 score matrix (B) and the substitution matrix (P).
* **[FastSP.jar](FastSP.jar):** Given to you. Please **DO NOT** modify or remove.
This file computes the SP-Score to measure accuracy of your alignment on part 2 (more details are given in the **Measurement of error** section).
* **[autocheck.py](autocheck.py):** Given to you. Please **DO NOT** modify. This file provides a test-suite for your program. After you complete the assignment, you will want to run `autocheck.py` to test your program by typing ```python3 autocheck.py```.

## Testing
As mentioned above, you can run ```python3 autocheck.py``` to automatically test your program.
There are 14 tests in total. The **first two tests** are sanity checks of the package distributed to you, so your program will always pass these first two tests. If it does not, please contact your TA to troubleshoot the problem.

The tests on ``autocheck.py`` will also be run by "Github Actions" everytime you push your code to the Github repository, so be sure to check the "Actions" tab of your repository to see if the tests passed.

Additionally, you can manually test your program using the test cases in [testing/test_data/checking](testing/test_data/checking) folder. 
To test part 1, you can use the `DP` folder. Run `threeway_align.py` with `-i` is one of the `.fas` file and `-g` is either -1 (for "g1") or -5 (for "g5"), then compare your output to the corresponding `.aln` file. Note that the correct way to compare the alignments is to compare their sum-pairwise scores. You can compute the alignment score by using the script `compute_pairwsie_score.py`. Type `compute_pairwise_score.py -h` to learn how to use it.

To test part 2, you can use the `highrate` and `lowrate` folders, each includes a set of `.fas` files, which are the *input sequences (unaligned)*, and the corresponding `.aln` files, which are the *true alignments*. Your goal is to match the output of your program to the true alignment as much as possible. For "lowrate" tests, the indel-rate is 0.01 and for "highrate" test, the indel-rate is 0.1.

For example, after you have implemented the `threeway_align` function, you would be able to call `threeway_align.py` as follow:

```python threeway_align.py -i testing/test_data/checking/lowrate/01.fas -o myOutput.aln -r 0.01```

Then you will need to compare `myOutput.aln` to `testing/test_data/checking/lowrate/01.aln` using `FastSP.jar` (details in the next section).

## Measurement of error: using [FastSP.jar](FastSP.jar) to compute the SP-Score
For part 1, your program must produce the optimal alignment that has the **highest sum of pairwise alignment scores**. 

For part 2, the gap penalty you computed will be supplied to the `threeway_align` function, and your program will be tested on synthetic data that was simulated using the BLOSUM62 model and the specified indel rate. 
Due to the stochastic nature of the simulated data, the alignment you produce is an "estimated alignment", and some level of uncertainties on the estimate is acceptable.

To evaluate the accuracy of your program, you will compare the output of your program to the true alignment. 
To do that, you will use **FastSP**<sup>2</sup> to compute the SP-Score. You are encouraged to read the paper to understand what SP-score is and how FastSP works.

The FastSP program, written in Java, is given to you in [FastSP.jar](FastSP.jar). You can run it to compute the SP-Score of your alignment (in addition to other accuracy measurements that will not be used here). If you have java installed, you can use the following command to run FastSP:

```java -jar FastSP.jar -r testing/test_data/checking/lowrate/01.aln -e myOutput.aln```

SP-score is a similarity measurement of two alignments. It takes values from 0 to 1, and the higher the SP-score, the more similar the two alignments are to each other. As mentioned before, you are not expected to produce perfect alignments. A good program will give alignments with high SP-scores, but not neccessarily obtain the perfect SP-Score of 1.0 for all tests. 

## Homework Deliverables
* **writeup.pdf:** A write-up that includes the following pieces of information:
    1. The recursive formula of the dynamic programming
    2. The base case(s) of the dynamic programming
    3. The asymptotic running time of the algorithm
    4. Your logic behind setting the gap penalty (i.e. part 2)
    * For each of the above items, you must derive or explain how you obtained the answer.

* **[threeway_align/todo.py](threeway_align/todo.py):** Your Python3 code for solving the problem
    * We have provided starter code, but you must complete the `threeway_align` and `compute_gap_penalty` functions (labeled with `TODO`).
    * Extra credits: improve the Dynamic Programming algorithm and add the improved version to `threeway_align_indel_rate` (labeled with `TODO: EXTRA CREDITS`).
    
## Grade Breakdown (100 Points)
* **Writeup: 40 Points**
    * *Correct recursive formula and base case(s):* 10 points
    * *Correct asymptotic running time:* 15 points
    * *Valid logic behind setting the gap penalty:* 15 points

* **Programming: 60 Points**
   * Part 1: 20 tests, 1.5 pts each. Time limit per test: 3 minutes
   * Part 2: 20 tests, 1.5 pts each. Time limit per test: 3 minutes. 
      * Low-rate: 10 tests, total 15 pts. Full credit: SP-Score >= 0.95. Half credit: SP-Score >= 0.8
      * High-rate: 10 tests, total 15 pts. Full credit: SP-Score >= 0.85. Half credit: SP-Score >= 0.8
      
* **Extra credits: up to 10 Points**      
   * Look beyond this dynamic programming algorithm to improve the accuracy of alignment or the running time of your program. You can come up with your own solution (be creative here!) or try a combination of the following ideas:
      * Implement *affine gap penalty* and figure out the penalties for "gap opening" and "gap extension" 
      * Come up with an *alternative objective function* that has better biological interpretation than the sum-of-pairwise function
      * Figure out a good way to *approximate* the optimal alignment score to save execution time
   * Note: extra credits will only be given to the solutions that lead to *substantial improvement* of either alignment accuracy or execution time (or both).
   
*Note:* the SP-Score thresholds are subject to change depending on the performance of the class, but if they change, they can only get *lower*, not *higher*.  

## References
1. Henikoff S, Henikoff JG. Amino acid substitution matrices from protein blocks. *Proceedings of the National Academy of Sciences*, 89(22):10915–10919, 1992. [doi:10.1073/pnas.89.22.10915](https://doi.org/10.1073/pnas.89.22.10915)
2. Mirarab S, Warnow T. FASTSP: linear time calculation of alignment accuracy. *Bioinformatics*, 27(23):3250–3258, 2011. [doi:10.1093/bioinformatics/btr553](https://doi.org/10.1093/bioinformatics/btr553)
