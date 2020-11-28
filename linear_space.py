# Followed https://en.wikipedia.org/wiki/Hirschberg's_algorithm for this one
import copy
import pdb
import numpy as np

LEFT = 1 
DIAG = 2
UP = 3

COL_MAP = {
	'A': 0,
	'C': 1,
	'T': 2,
	'G': 3,
	'-': 4
}

def get_max_score(scoreL, scoreR):
    max_index = 0
    max_sum = float('-Inf')
    for i, (l, r) in enumerate(zip(scoreL, scoreR[::-1])):
        # calculate the diagonal maximum index
        if sum([l, r]) > max_sum:
            max_sum = sum([l, r])
            max_index = i
    return max_index 



def NWScore(self, seqa, seqb):

    lena = len(seqa)
    lenb = len(seqb)
    pre_row = [0] * (lenb + 1)
    cur_row = [0] * (lenb + 1)

    for j in range(1, lenb + 1):
        pre_row[j] = pre_row[j - 1] + self.score_matrix.iloc[4]['A']

    for i in range(1, lena + 1):
        cur_row[0] = self.score_matrix.iloc[0]['-'] + pre_row[0]
        for j in range(1, lenb + 1):
            cur_row[j] = max(pre_row[j - 1] + self.score_matrix.iloc[COL_MAP[seqa[i-1]]][seqb[j-1]], 
                            pre_row[j] + self.score_matrix.iloc[0]['-'], 
                            cur_row[j - 1] + self.score_matrix.iloc[4]['A'])

        pre_row = cur_row
        cur_row = [0] * (lenb + 1)

    return pre_row


def dynamicProgramming(self, x, y):
    # M records is the score array
    # Path stores the path information, inside of Path:
    M = np.zeros((len(x) + 1, len(y) + 1))
    Path = np.empty((len(x) + 1, len(y) + 1), dtype=object)

    for i in range(1, len(y) + 1):
        M[0][i] = M[0][i-1] + self.score_matrix.iloc[4]['A']
        Path[0][i] = LEFT
    for j in range(1, len(x) + 1):
        M[j][0] = M[j-1][0] + self.score_matrix.iloc[0]['-']
        Path[j][0] = UP

    for i in range(1, len(x) + 1):
        for j in range(1, len(y) + 1):
        
            M[i][j] = max(M[i-1][j-1] + self.score_matrix.iloc[COL_MAP[x[i-1]]][y[j-1]], M[i-1][j] + self.score_matrix.iloc[4]['A'], M[i][j-1] + self.score_matrix.iloc[0]['-'])
            if M[i][j] == M[i-1][j-1] + self.score_matrix.iloc[COL_MAP[x[i-1]]][y[j-1]]:
                Path[i][j] =  DIAG
            elif M[i][j] == M[i-1][j] + self.score_matrix.iloc[4]['A']:
                Path[i][j] = UP
            else:
                Path[i][j] = LEFT

    row = []
    column= []
    i = len(x)
    j = len(y)
    while Path[i][j]:
        if Path[i][j] == DIAG:
            row.insert(0, y[j-1])
            column.insert(0, x[i-1])

            i -= 1
            j -= 1
        elif Path[i][j] == UP:
            row.insert(0, '-')
            column.insert(0, x[i-1])
            i -= 1
        elif Path[i][j] == LEFT:
            column.insert(0, '-')
            row.insert(0, y[j-1])
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

        xmid = len(x)// 2
        ylen = len(y)

        scoreL = NWScore(self, x[:xmid], y)
        scoreR = NWScore(self, x[xmid:][::-1], y[::-1])
        ymid = get_max_score(scoreL, scoreR)

        row_l, column_u = Hirschberg(self, x[:xmid], y[:ymid])
        row_r, column_d = Hirschberg(self, x[xmid:], y[ymid:])
        row = row_l + row_r
        column = column_u + column_d 
    return row, column


class LinearSpaceAlignment():

    def __init__(self, seq1, seq2, score_matrix, mode):
        self.seq1 = seq1
        self.seq2 = seq2
        self.mode = mode
        self.row_size = len(seq1) + 1
        self.col_size = len(seq2) + 1
        self.score_matrix = score_matrix
        z,w = Hirschberg(self, self.seq1, self.seq2)
        print(w)
        print(z)




if __name__ == "__main__":
    LinearSpaceAlignment()