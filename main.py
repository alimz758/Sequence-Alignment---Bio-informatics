import sys
import pandas as pd
import time
from global_local_alignment import SequenceInit
from linear_space import LinearSpaceAlignment
#from test_linear import LinearSpaceAlignment


ACCEPTED_MODES = ["global", "local"]
ACCEPTED_CHARS = ['A', 'C', 'T', 'G']


def get_score_matrix():
    file = sys.argv[4]
    if not file.endswith('.csv'):
        raise ValueError("Only .csv is accepted for score_matrix")
    return pd.read_csv(file)

def get_mode():
    mode = sys.argv[1]
    if not mode in ACCEPTED_MODES:
        raise ValueError("Only " + str(ACCEPTED_MODES) + " are accepted")
    return mode

    
def get_sequence():
    
    file1 = sys.argv[2]
    file2 = sys.argv[3]
    if not file1.endswith('.txt') and not file2.endswith('.txt'):
        raise ValueError("Only .txt is accepted for sequence files")
    
    seq1 = open(file1, "r").read()
    seq2 = open(file2, "r").read()
    if any(c not in ACCEPTED_CHARS for c in seq1) or any(c not in ACCEPTED_CHARS for c in seq2): 
        raise ValueError("Only " + str(ACCEPTED_CHARS) + " characters are accepted")
    
    return seq1, seq2


def get_alignments():
    mode = get_mode()
    seq1, seq2 = get_sequence()
    score_matrix_df = get_score_matrix()
    print("\n================== Results ==================")
    print("Sequence1: ", seq1)
    print("Sequence2: ", seq2)
    print("Total sequences length is: ", len(seq1) + len(seq2))
   
     # Needleman–Wunsch_algorithm
    t0 = time.time()
    sequence = SequenceInit(seq1, seq2, score_matrix_df, mode)
    t1 = time.time()
    sequence.print_results(t1, t0)
    # Hirschberg
    ls = LinearSpaceAlignment(seq1, seq2, score_matrix_df, mode)
    ls.print_results()
   
    
def main():
    get_alignments()


if __name__ == "__main__":
    main()