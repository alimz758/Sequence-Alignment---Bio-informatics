# Sequence-Alignment---Bio-informatics

Course project for CS121 at UCLA.

For a full description please take a look at `specs.pdf`, project#2.

This program can run for gloabl alignment and local alignment on two different alogirthms. `Needleman–Wunsch_algorithm` and `Hirschberg` for space efficiency

## How to Run

To run this programm you need to specify 4 things.

1. Mode 

2. Sequence1 as a text file

3. Sequence2 as a text file

4. Score matrix as a csv file

### Modes

There are 4 different modes to get sequence alignments:

1. `global`: Get all global alignments using `Needleman–Wunsch_algorithm`

2. `local`: Get all local alignments using `Needleman–Wunsch_algorithm`

3. `middle-global`: Get global alignment using `Hirschberg`

4. `middle-local`: Get local alignment using `Hirschberg`

For instance this would return all global alignments using `Needleman–Wunsch_algorithm`:

```
> python main.py global seq1.txt seq2.txt score_matrix.csv
Total sequence length is:  7
[[( 0, 0) (-1, 'L') (-2, 'L') (-3, 'L') (-4, 'L')]
 [(-1, 'U') ( 2, 'D') ( 1, 'L') ( 0, 'L') (-1, 'LD')]
 [(-2, 'U') ( 1, 'U') ( 4, 'D') ( 3, 'D') ( 2, 'L')]
 [(-3, 'U') ( 0, 'DU') ( 3, 'U') ( 3, 'D') ( 5, 'D')]]
Max global alignment Score is: 5
Alignmnet # 1
ATTA
AT-A

Alignmnet # 2
ATTA
A-TA
```


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