import numpy as np
import math


SEQ1 = 'ACG' #Sequence 1 (Side Sequence)
SEQ2 = 'AACT' #Sequence 2 (Top Sequence)
SCORE_MATRIX = {
'A': {'A': 1, 'C':0, 'G':0, 'T':0}, 
'C': {'A':0, 'C': 1, 'G':0, 'T':0}, 
'G': {'A':0, 'C':0, 'G': 1, 'T':0},
'T': {'A':0, 'C':0, 'G':0, 'T': 1},
'-': -1
}
#back pointers
LEFT = 1 
DIAG = 2
UP = 3
SOURCE = 4


def set_first_row_col_to_zero(self, isglobal):
    #Set the last col of every row
    for i in range(self.row_size):
        if isglobal:
            self.dp_matrix[i][0] = [self.indel*i, [LEFT]] 
        else: 
            self.dp_matrix[i][0] = [0, [SOURCE]] 
            
    #Set the last row of every col    
    for j in range(self.col_size):
        if isglobal:
            self.dp_matrix[0][j] = [self.indel*j, [UP]]
        else: 
            self.dp_matrix[0][j] = [0, [SOURCE]]               

def discover_global_paths(self, i, j, dest_x=0, dest_y=0, path=''):
    if i == dest_x and j==dest_y: 
        self.alignment_paths.append(path) 
        return 

    current_dir_options = len(self.dp_matrix[i][j][1]) 
    while current_dir_options<=1: 
        current_dir = self.dp_matrix[i][j][1][0]  
        path = path + str(current_dir) 
        if current_dir == LEFT: 
            j=j-1
        elif current_dir == DIAG:
            i=i-1
            j=j-1
        elif current_dir == UP:
            i=i-1
        current_dir_options = len(self.dp_matrix[i][j][1])
        if i == dest_x and j==dest_y or current_dir == SOURCE:
            self.alignment_paths.append(path)
            return 
        
    if current_dir_options>1:
        for dir_option in range(current_dir_options):
            current_dir = self.dp_matrix[i][j][1][dir_option] 
            tmp_path = path + str(current_dir)
            if current_dir == LEFT:
                n_i = i
                n_j=j-1
            elif current_dir == DIAG:
                n_i=i-1
                n_j=j-1
            elif current_dir == UP:
                n_i=i-1
                n_j = j
            discover_global_paths(self, n_i, n_j, dest_x, dest_y, tmp_path)
            
def discover_local_paths(self, i, j, dest_x=0, dest_y=0, path=''):
    if i == dest_x and j==dest_y: 
        self.alignment_paths.append(path) 
        return 

    current_dir_options = len(self.dp_matrix[i][j][1]) 
    while current_dir_options<=1: 
        current_dir = self.dp_matrix[i][j][1][0] 
        path = path + str(current_dir) 
        if current_dir == LEFT: 
            j=j-1
        elif current_dir == DIAG:
            i=i-1
            j=j-1
        elif current_dir == UP:
            i=i-1
        current_dir_options = len(self.dp_matrix[i][j][1]) 
        if i == dest_x and j==dest_y or current_dir == SOURCE:
            self.alignment_paths.append(path)
            return 
        
    if current_dir_options>1:
        for dir_option in range(current_dir_options):
            current_dir = self.dp_matrix[i][j][1][dir_option]
            tmp_path = path + str(current_dir)
            if current_dir == LEFT:
                n_i = i
                n_j=j-1
            elif current_dir == DIAG:
                n_i=i-1
                n_j=j-1
            elif current_dir == UP:
                n_i=i-1
                n_j = j
            elif current_dir == SOURCE:
                self.alignment_paths.append(path)
                return
            discover_global_paths(self, n_i, n_j, dest_x, dest_y, tmp_path)
             

def get_all_alignments(self, max_score_i, max_score_j):
    aln_count = 0
    for elem in self.alignment_paths:
        i = max_score_i-1
        j = max_score_j-1
        side_aln = ''
        top_aln = ''
        step = 0
        aln_info = []
        for n_dir_c in range(len(elem)):
            n_dir = elem[n_dir_c]
            score = self.dp_matrix[i+1][j+1][0]
            step = step + 1
            aln_info.append([step,score,n_dir])
            if n_dir == '2':
                side_aln = side_aln +  self.seq1[i]
                top_aln = top_aln +  self.seq2[j]
                i=i-1
                j=j-1
            elif n_dir == '1':
                side_aln = side_aln + '-'
                top_aln = top_aln +  self.seq2[j]
                j=j-1
            elif n_dir == '3':
                side_aln = side_aln +  self.seq1[i]
                top_aln = top_aln + '-'
                i=i-1
        aln_count = aln_count + 1
        self.all_alignments.append([top_aln[::-1],side_aln[::-1],elem,aln_info, aln_count])

def print_alignments(self, type):
    print('Total Number of {} Alignments: {}'.format(type, len(self.all_alignments)))
    print('Max Score: '+str(self.all_alignments[0][3][0][1])+'\n')
    print("Top string {}, bottom string: {}".format(self.seq2, self.seq1))
    for alignment in self.all_alignments:
        print(alignment[0]+ '\n' + alignment[1]+ '\n')
        
class SequenceAlignment(object):
    
    def __init__(self, seq1, seq2, score_matrix):
        self.seq1 = seq1
        self.seq2 = seq2
        self.row_size = len(seq1) + 1
        self.col_size = len(seq2) + 1
        self.score_matrix = score_matrix
        self.indel = self.score_matrix['-']
        self.dp_matrix = [[[[None] for i in range(2)] for i in range( self.col_size)] for i in range(self.row_size)] # create dp_matrix to keep rack of score and directions [SCORE, [dir1, dir2]]
        self.all_alignments = []
        self.alignment_paths = []
      
    #Create the DP matrix score and direction for the global alignment  
    def create_global_alignment(self):
        set_first_row_col_to_zero(self, True)
        for i in (range(1, self.row_size)):
            for j in range(1, self.col_size):
                char1 = self.seq1[i-1]
                char2 = self.seq2[j-1]
                match_or_mismatch_score = self.score_matrix[char1][char2]
                score_list = [self.dp_matrix[i][j-1][0] + self.indel , self.dp_matrix[i-1][j-1][0] + match_or_mismatch_score, self.dp_matrix[i-1][j][0] + self.indel]
                max_score = max(score_list)
                self.dp_matrix[i][j] = [max_score, [i+1 for i,v in enumerate(score_list) if v==max(score_list)]]
        #max_score = self.dp_matrix[i][j][0]
        max_score_i = i
        max_score_j = j
        discover_global_paths(self, max_score_i, max_score_j)
        get_all_alignments(self, max_score_i, max_score_j)
        print_alignments(self, "Global")
        
    #Create the DP matrix score and direction for the local alignment  
    def create_local_alignment(self):
        set_first_row_col_to_zero(self, False)
        best_score = -math.inf
        for i in (range(1, self.row_size)):
            for j in range(1, self.col_size):
                char1 = self.seq1[i-1]
                char2 = self.seq2[j-1]
                match_or_mismatch_score = self.score_matrix[char1][char2]
                score_list = [self.dp_matrix[i][j-1][0] + self.indel , self.dp_matrix[i-1][j-1][0] + match_or_mismatch_score, self.dp_matrix[i-1][j][0] + self.indel, 0]
                max_score = max(score_list)
                self.dp_matrix[i][j] = [max_score, [i+1 for i,v in enumerate(score_list) if v==max(score_list)]]
                if max_score > best_score:
                    best_score = max_score
                    max_score_i = i
                    max_score_j = j
        #max_score = self.dp_matrix[i][j][0]
        #print(self.dp_matrix)
        discover_local_paths(self, max_score_i, max_score_j)
        get_all_alignments(self, max_score_i, max_score_j)
        print_alignments(self, "Local")
    
        
                
                    
            
def main():
   sq = SequenceAlignment(SEQ1, SEQ2, SCORE_MATRIX)
   # Create a global alignment
   #sq.create_global_alignment()
   sq.create_local_alignment()


if __name__ == "__main__":
    main()