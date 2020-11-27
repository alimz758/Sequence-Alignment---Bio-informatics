import sys
import numpy as np

COL_MAP = {
	'A': 0,
	'C': 1,
	'T': 2,
	'G': 3,
	'-': 4
}


def get_optimal_point(self, i, j):
	max_score = float('-inf')
	for r in range(1,self.row):
		for c in range(1, self.column):
			if self.path_matrix[r][c][0]>=max_score:
				i = r
				j = c
				max_score = self.path_matrix[r][c][0]
	return i, j, max_score


def get_max_score(self):

	return self.path_matrix[self.row-1][self.column-1][0]
 
			
def get_sequences(self, i, j, aligned_seq1 = "", aligned_seq2 = ""):
	
	if self.path_matrix[i][j][1]==0:
		print(aligned_seq1)
		print(aligned_seq2)
		print()
		return 

	if len(self.path_matrix[i][j][1])>1:
		paths = self.path_matrix[i][j][1]
		for n in range(len(paths)):
			self.path_matrix[i][j][1] = paths[n]
			get_sequences(self, i, j, aligned_seq1, aligned_seq2)

	else:	
		if self.path_matrix[i][j][1] == 'D':
			aligned_seq1 = self.seq1[j-1] + aligned_seq1
			aligned_seq2 = self.seq2[i-1] + aligned_seq2
			i = i-1
			j = j-1

		elif self.path_matrix[i][j][1]=='U':
			aligned_seq1 = "-" + aligned_seq2
			aligned_seq2 = self.seq2[i-1] + aligned_seq2
			i = i-1

		else:
			aligned_seq2 = "-" + aligned_seq2
			aligned_seq1 = self.seq1[j-1] + aligned_seq1
			j = j-1

		get_sequences(self, i, j, aligned_seq1, aligned_seq2)			


def get_path(self, i, j, mode):

	paths = ''
	match_mismatch_score = self.score.iloc[COL_MAP[self.seq2[i-1]]][self.seq1[j-1]]
	diag = self.path_matrix[i-1][j-1][0] + match_mismatch_score
	up = self.path_matrix[i-1][j][0] + self.score.iloc[4]['A']
	left = self.path_matrix[i][j-1][0] + self.score.iloc[0]['-']

	self.path_matrix[i][j][0] = max(diag, up, left)

	if self.path_matrix[i][j][0]==left:
		paths = paths + 'L'

	if self.path_matrix[i][j][0]==diag:
		paths = paths + 'D'

	if self.path_matrix[i][j][0]==up:
		paths = paths + 'U'

	self.path_matrix[i][j][1] = paths

	if i==self.row-1 and j==self.column-1:
		print(self.path_matrix)
		score = 0
		if mode == "global":
			score = get_max_score(self)
		elif mode == "local":
			i, j, score = get_optimal_point(self, i, j)

		print('{} alignment Score is: {}'.format(mode, score))
		get_sequences(self, i, j)
		
	elif j<self.column-1:
		get_path(self, i ,j+1, mode)
		
	else:
		get_path(self, i+1 ,1, mode)


def init_path_matrix(self):
	
	matrix = np.zeros([self.row, self.column], dtype='i,O') 
	
	for i in range(1,self.column):
		matrix[0][i][0] = i * self.score.iloc[0]['-']
		matrix[0][i][1] = 'L'

	for i in range(1,self.row):
		matrix[i][0][0] = i * self.score.iloc[4]['A']
		matrix[i][0][1] = 'U'
		
	return matrix


class SequenceInit(object):
	
	def __init__(self, seq1, seq2, score, mode):
		self.seq1 = seq1
		self.seq2 = seq2
		self.score = score
		self.column = len(self.seq1)+1
		self.row = len(self.seq2)+1
		self.path_matrix = init_path_matrix(self)
		get_path(self, 1, 1, mode)	