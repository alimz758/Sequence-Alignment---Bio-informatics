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

> python main.py global seq1.txt seq2.txt score_matrix.csv