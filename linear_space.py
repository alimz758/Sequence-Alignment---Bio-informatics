# Followed https://en.wikipedia.org/wiki/Hirschberg's_algorithm for this one
import copy
import pdb
import numpy as np
from sequenence_init import SequenceInit

LEFT = 1 
DIAG = 2
UP = 3
SOURCE = 4

def get_max_score(scoreL, scoreR):
    max_index = 0
    max_sum = float('-Inf')
    for i, (l, r) in enumerate(zip(scoreL, scoreR[::-1])):
        # calculate the diagonal maximum index
        if sum([l, r]) > max_sum:
            max_sum = sum([l, r])
            max_index = i
    return max_index 

def test(self, x, y):
    M = np.zeros((len(x)+1, len(y)+1))
    for i in range(len(x)+1):
        M[i][0] = i * self.indel
    for j in range(len(y)+1 ):
        M[0][j] = j * self.indel
    for i in range(1, len(x)+1):
        for j in range(1, len(y)+1):
            M[i][j] = max(M[i-1][j-1] + self.score_matrix[x[i-1]][y[i-1]], M[i-1][j] + self.indel, M[i][j-1] + self.indel)
    
    AlignmentA = ''
    AlignmentB = ''
    i = len(x)
    j = len(y)
    while i > 0 or j > 0:
        if i > 0 and j > 0 and (M[i][j] == M[i-1][j-1] + self.score_matrix[x[i-1]][y[j-1]]):
            AlignmentA = x[i-1] + AlignmentA
            AlignmentB = y[j-1] + AlignmentB
            i = i - 1
            j = j - 1
        elif i > 0 and M[i][j] == M[i-1][j] + self.indel:
            AlignmentA = x[i-1] + AlignmentA
            AlignmentB = '-' + AlignmentB
            i = i - 1
        else:
            AlignmentA = '-' + AlignmentA
            AlignmentB = y[j-1] + AlignmentB
            j = j - 1
    print(AlignmentA)
    return AlignmentA, AlignmentB 


def NWScore(self, x, y):

    row = y
    column = x
    minLen = len(y)
    prev = [0 for i in range(minLen + 1)]
    current = [0 for i in range(minLen + 1)]
    for i in range(1, minLen + 1):
        prev[i] = prev[i-1] + self.indel
    
    current[0] = 0
    for j in range(1, len(column) + 1):
        current[0] += self.indel
        for i in range(1, minLen + 1):
            if row[i-1] == column[j-1]:
                try:
                    current[i] = max(current[i-1] + self.indel, prev[i-1] + self.score_matrix[x[j-1]][y[i-1]], prev[i] + self.indel)
                except:
                    pdb.set_trace()
            else:
                current[i] = max(current[i-1] + self.indel, prev[i-1] + self.score_matrix[x[j-1]][y[i-1]], prev[i] + self.indel)
        prev = copy.deepcopy(current) 

    return current 


        
def dynamicProgramming(self, x, y):
    # M records is the score array
    # Path stores the path information, inside of Path:
    M = np.zeros((len(x) + 1, len(y) + 1))
    Path = np.empty((len(x) + 1, len(y) + 1), dtype=object)

    for i in range(1, len(y) + 1):
        M[0][i] = M[0][i-1] + self.indel
        Path[0][i] = LEFT
    for j in range(1, len(x) + 1):
        M[j][0] = M[j-1][0] + self.indel
        Path[j][0] = UP
    
    for i in range(1, len(x) + 1):
        for j in range(1, len(y) + 1):
            if x[i-1] == y[j-1]:
                M[i][j] = max(M[i-1][j-1] + self.score_matrix[x[i-1]][y[j-1]], M[i-1][j] + self.indel, M[i][j-1] + self.indel)
                if M[i][j] == M[i-1][j-1] + self.score_matrix[x[i-1]][y[j-1]]:
                    Path[i][j] =  DIAG
                elif M[i][j] == M[i-1][j] + self.indel:
                    Path[i][j] = UP
                else:
                    Path[i][j] = LEFT
            else:
                M[i][j] = max(M[i-1][j-1] + self.score_matrix[x[i-1]][y[j-1]], M[i-1][j] + self.indel, M[i][j-1] + self.indel)
                if M[i][j] == M[i-1][j-1] + self.score_matrix[x[i-1]][y[j-1]]:
                    Path[i][j] =  DIAG
                elif M[i][j] == M[i-1][j] + self.indel:
                    Path[i][j] = UP
                else:
                    Path[i][j] = LEFT

    row = []
    column= []
    middle = []
    i = len(x)
    j = len(y)
    while Path[i][j]:
        if Path[i][j] == DIAG:
            row.insert(0, y[j-1])
            column.insert(0, x[i-1])
            if x[i-1] == y[j-1]:
                middle.insert(0, '|')
            else:
                middle.insert(0, ':')
            i -= 1
            j -= 1
        elif Path[i][j] == UP:
            row.insert(0, '-')
            column.insert(0, x[i-1])
            middle.insert(0, 'x')
            i -= 1
        elif Path[i][j] == LEFT:
            column.insert(0, '-')
            row.insert(0, y[j-1])
            middle.insert(0, 'x')
            j -= 1
            
    return row, column 

def Hirschberg(self, x, y):
    row = ""
    column = ""
    if len(x) == 0:
        column = '-' * len(y)
        row = y
    elif len(y) == 0:
        column += x
        row += '-' * len(x)
    elif len(x) == 1 or len(y) == 1:
        row,column = dynamicProgramming(self, x, y)
        row, column = map(lambda x: "".join(x), [row, column]) 
    else:
        
        xmid = len(x)/ 2
        ylen = len(y)

        scoreL = NWScore(self, x[:xmid], y)
        scoreR = NWScore(self, x[xmid:][::-1], y[::-1])
        ymid = get_max_score(scoreL, scoreR)
        
        row_l, column_u = Hirschberg(self, x[:xmid], y[:ymid])
        row_r, column_d = Hirschberg(self, x[xmid:], y[ymid:])
        row = row_l + row_r
        column = column_u + column_d 
    return row, column


class LinearSpaceAlignment(SequenceInit):
    
    def __init__(self, seq1, seq2, score):
        super(LinearSpaceAlignment, self).__init__(seq1, seq2, score)
        z,w = Hirschberg(self, self.seq1, self.seq2)
        print(z)
        print(w)
 
