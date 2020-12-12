# Sequence-Alignment---Bio-informatics

Course project for CS121 at UCLA.

For a full description please take a look at `specs.pdf`, project#2.

This program can run for global alignment and local alignments utilizing two different alogirthms: `Needleman–Wunsch_algorithm` and `Hirschberg` for space efficiency

## How to Run

> python main.py  \<MODE> <SEQ1.txt> <SEQ2.txt> <SCORE_MATRIX.csv> 

To run this programm you need to specify 4 things in the above order.

1. Mode 

2. Sequence1 as a text file

3. Sequence2 as a text file

4. Score matrix as a csv file

### Modes

There are 2 different modes to get sequence alignments:

1. `global`: Get all global alignments using `Needleman–Wunsch` and `Hirschberg`

2. `local`: Get all local alignments using `Needleman–Wunsch` and `Hirschberg`

### SEQ1.txt

Text file containing the Sequence1 characters(A, C, T, G are acceptable) 

### SEQ2.txt

Text file containing the Sequence2 characters(A, C, T, G are acceptable) 

### Score Matrix

The score matrix is a `csv` file that follows the same format as desired by the project. If you wish to change the scores, you can directly modify the csv file.

You can see its format below:

```
 ,  A,  C,  T,  G,  -
A,  2,  -1, -1, -1, -1
C,  -1, 2,  -1, -1, -1
T,  -1, -1, 2,  -1, -1
G,  -1, -1, -1, 2,  -1
-,  -1, -1, -1, -1

```

For instance this would return all local alignments using both algorithms:

> python main.py local seq1.txt seq2.txt score_matrix.csv
```
================== Results ==================
Sequence1:  ATTCG
Sequence2:  ATG
Total sequences length is:  8


==========   Needleman–Wunsch  ==========
[[(0, 0) (0, 0) (0, 0) (0, 0) (0, 0) (0, 0)]
 [(0, 0) (2, 'D') (0, 'L') (0, '') (0, '') (0, '')]
 [(0, 0) (0, 'U') (4, 'D') (2, 'LD') (0, 'L') (0, '')]
 [(0, 0) (0, '') (2, 'U') (3, 'D') (1, 'LD') (2, 'D')]]
Max local alignment Score is: 4
Total run time in seconds:  0.016
Alignmnet # 1
AT
AT

Memory Usage(using `tracemalloc`, the first stats):
 size=7930 B (+7930 B), count=109 (+109), average=73 B



==========   Hirschberg  ==========
Max local alignment Score is: 4
Total run time in seconds:  0.028
Alignmnet # 1
AT
AT

Memory Usage(using `tracemalloc`, the first stats):
 size=1408 B (+1408 B), count=3 (+3), average=469 B
```
