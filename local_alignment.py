#!/usr/bin/python
import sys
import numpy as np


def get_max_score(self,seq1, seq2):
    if (len(seq1)== 0):
        return
    score = 0
    for i in range(len(seq1)):
        if seq1[i] == seq2[i]:
            score = score + self.match
        else:
            if seq1[i]=='-' or seq2[i]=='-':
                score = score + self.d
            else:
                score = score + self.mismatch

    print ('Local Alignment Score is: ', score)

def get_sequences(self, F, i, j, alignmented_seq1 = "", alignmented_seq2 = ""):

	if F[i][j][0]==0:
		get_max_score(self, alignmented_seq1, alignmented_seq2)	
		print (alignmented_seq1)
		print (alignmented_seq2)	
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

		get_sequences(self, F, i, j, alignmented_seq1, alignmented_seq2)			


def local_alignment(self, F, i, j):

	value = self.match
	if self.seq2[i-1] != self.seq1[j-1]:
		value = self.mismatch  

	directions = ''	

	diag = F[i-1][j-1][0] + value
	up = F[i-1][j][0] + self.d
	left = F[i][j-1][0] + self.d

	F[i][j][0]	= max(diag, up, left, 0)
	
	if F[i][j][0]==diag:
		directions = directions + 'D'

	if F[i][j][0]==up:	
		directions = directions + 'U'
	
	if F[i][j][0]==left:
		directions = directions + 'L'	

	F[i][j][1] = directions
		
	if i==self.row-1 and j==self.column-1:

		major = float('inf')
		s_i = 0
		s_j = 0 
		for r in range(1,self.row):
			for c in range(1, self.column):
				if F[r][c][0]>=major:
					s_i = i
					s_j = j
					i = r
					j = c
					major = F[r][c][0]

		get_sequences(self, F, i, j)
		get_sequences(self, F, s_i, s_j)
		return

	if j<self.column-1:
		local_alignment(self, F, i ,j+1)	
	else:	
		local_alignment(self, F, i+1 ,1)


class LocalAlignment():
    
    def __init__(self, seq1, seq2, score):

        self.seq1 = seq1
        self.seq2 = seq2
        self.match = 2
        self.mismatch = -1	
        self.d = -1
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
        local_alignment(self, F, 1, 1)	