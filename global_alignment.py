#!/usr/bin/python
import sys
import numpy as np


def get_max_score(self,seq1, seq2):
    score = 0
    for i in range(len(seq1)):
        if seq1[i] == seq2[i]:
            score = score + self.match
        else:
            if seq1[i]=='-' or seq2[i]=='-':
                score = score + self.d
            else:
                score = score + self.mismatch

    print ('Global Alignment Score is: ', score)
 
			
def get_sequences(self,F, i, j, alignmented_seq1 = "", alignmented_seq2 = ""):
    if F[i][j][1]==0:
        get_max_score(self,alignmented_seq1, alignmented_seq2)
        print (alignmented_seq1)
        print (alignmented_seq2)
        print()
        return

    if len(F[i][j][1])>1:
        directions = F[i][j][1]
        for n in range(len(directions)):
            F[i][j][1] = directions[n]
            get_sequences(self, F, i, j, alignmented_seq1, alignmented_seq2)

    else:	
        if F[i][j][1] == 'D':
            alignmented_seq1 = self.seq1[j-1] + alignmented_seq1
            alignmented_seq2 = self.seq2[i-1] + alignmented_seq2
            i = i-1
            j = j-1

        elif F[i][j][1]=='U':
            alignmented_seq1 = "-" + alignmented_seq2
            alignmented_seq2 = self.seq2[i-1] + alignmented_seq2
            i = i-1

        else:
            alignmented_seq2 = "-" + alignmented_seq2
            alignmented_seq1 = self.seq1[j-1] + alignmented_seq1
            j = j-1

        get_sequences(self,F, i, j, alignmented_seq1, alignmented_seq2)			


def global_alignment(self,F, i, j):

    directions = ''

    value = self.match
    if self.seq2[i-1] != self.seq1[j-1]:
        value = self.mismatch  

    diag = F[i-1][j-1][0] + value
    up = F[i-1][j][0] + self.d
    left = F[i][j-1][0] + self.d

    F[i][j][0] = max(diag, up, left)

    if F[i][j][0]==left:
        directions = directions + 'L'

    if F[i][j][0]==diag:
        directions = directions + 'D'

    if F[i][j][0]==up:
        directions = directions + 'U'

    F[i][j][1] = directions

    if i==self.row-1 and j==self.column-1:
        get_sequences(self,F, i, j)
        return

    elif j<self.column-1:
        global_alignment(self,F, i ,j+1)
    else:
        global_alignment(self,F, i+1 ,1)



class GlobalAlignment():

    def __init__(self, seq1, seq2, score):
        self.seq1 = seq1
        self.seq2 = seq2
        self.match = 2
        self.mismatch = -1	
        self.d = -2
        self.score = score
        # zeros column and row at the beginning of matrix F
        self.column = len(self.seq1)+1
        self.row = len(self.seq2)+1

        F = np.zeros([self.row, self.column], dtype='i,O') 
        # Adding zeros 
        for i in range(1,self.column):
            F[0][i][0] = i * self.d
            F[0][i][1] = 'L'

        for i in range(1,self.row):
            F[i][0][0] = i * self.d
            F[i][0][1] = 'U'
        global_alignment(self, F, 1, 1)	