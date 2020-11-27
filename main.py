import sys
from global_alignment import GlobalAlignment
from local_alignment import LocalAlignment
from linear_space import LinearSpaceAlignment

SEQ1 = "ATTA" 
SEQ2 = "ATA"
SCORE_MATRIX = {
    'A': {'A': 2, 'C':-1, 'G':-1, 'T':-1}, 
    'C': {'A':-1, 'C': 2, 'G':-1, 'T':-1}, 
    'G': {'A':-1, 'C':-1, 'G': 2, 'T':-1},
    'T': {'A':-1, 'C':-1, 'G':-1, 'T': 2},
    '-': -2
}
ACCEPTED_MODES = ["global", "local", "middle"]

# def read_from_file():
#     f = open("sequences_info.txt", "r")
#     info = f.read()
#     print(info(1))

def get_mode():
    mode = sys.argv[1]
    if not mode in ACCEPTED_MODES:
        raise ValueError("Only " + str(ACCEPTED_MODES) + " are accepted")
    return mode

def get_alignments(mode):

    if mode == "global":
        get_global()
    elif mode == "local":
        get_local()
    elif mode == "middle":
        get_middle()
           
def get_global():
    GlobalAlignment(SEQ1, SEQ2, SCORE_MATRIX)

def get_local():
    LocalAlignment(SEQ1, SEQ2, SCORE_MATRIX)

def get_middle():
    LinearSpaceAlignment(SEQ1, SEQ2, SCORE_MATRIX)
    
    
# def get_local():
    
        
def main():
    # read_from_file()
    mode = get_mode()
    get_alignments(mode)


if __name__ == "__main__":
    main()