import sys
import numpy as np
import time
import pandas as pd
import tracemalloc


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
			if self.path_matrix[r][c][0] >= max_score:
				i = r
				j = c
				max_score = self.path_matrix[r][c][0]
	return i, j, max_score


def get_max_score(self):

	return self.path_matrix[self.row-1][self.column-1][0]
 
			
def get_sequences(self, i, j, aligned_seq1 = "", aligned_seq2 = ""):
	
	if (self.mode == "global" and self.path_matrix[i][j][1]==0) or (self.mode == "local" and self.path_matrix[i][j][0]==0):
		self.alignments.append([aligned_seq1, aligned_seq2])
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

		elif self.path_matrix[i][j][1] == 'U':
			aligned_seq1 = "-" + aligned_seq1
			aligned_seq2 = self.seq2[i-1] + aligned_seq2
			i = i-1

		else:
			aligned_seq2 = "-" + aligned_seq2
			aligned_seq1 = self.seq1[j-1] + aligned_seq1
			j = j-1

		get_sequences(self, i, j, aligned_seq1, aligned_seq2)			


def get_path(self, i, j):

	paths = ''
	match_mismatch_score = self.score_matrix.iloc[COL_MAP[self.seq2[i-1]]][self.seq1[j-1]]
	diag = self.path_matrix[i-1][j-1][0] + match_mismatch_score
	up = self.path_matrix[i-1][j][0] + self.score_matrix.iloc[4]['A']
	left = self.path_matrix[i][j-1][0] + self.score_matrix.iloc[0]['-']

	if self.mode == "global":
		self.path_matrix[i][j][0] = max(diag, up, left)
	else:
		self.path_matrix[i][j][0] = max(diag, up, left, 0)

	if self.path_matrix[i][j][0]==left:
		paths = paths + 'L'

	if self.path_matrix[i][j][0]==diag:
		paths = paths + 'D'

	if self.path_matrix[i][j][0]==up:
		paths = paths + 'U'

	self.path_matrix[i][j][1] = paths

	if i==self.row-1 and j==self.column-1:
		if self.mode == "global":
			self.max_score = get_max_score(self)
		elif self.mode == "local":
			i, j, self.max_score = get_optimal_point(self, i, j)
		get_sequences(self, i, j)
		
	elif j<self.column-1:
		get_path(self, i ,j+1)

	else:
		get_path(self, i+1 ,1)


def init_path_matrix(self):
	
	matrix = np.zeros([self.row, self.column], dtype='i,O') 
	if self.mode == "local":
		return matrix
	for i in range(1, self.column):
		matrix[0][i][0] = i * self.score_matrix.iloc[0]['-']
		matrix[0][i][1] = 'L'

	for i in range(1, self.row):
		matrix[i][0][0] = i * self.score_matrix.iloc[4]['A']
		matrix[i][0][1] = 'U'
		
	return matrix


class SequenceInit(object):
	
	def __init__(self, seq1, seq2, score_matrix, mode):
	 
		self.seq1 = seq1
		self.seq2 = seq2
		self.mode = mode
		self.score_matrix = score_matrix
		self.column = len(self.seq1)+1
		self.row = len(self.seq2)+1
		self.max_score = 0
		self.alignments = []
		self.path_matrix = init_path_matrix(self)
		tracemalloc.start()
		snapshot1 = tracemalloc.take_snapshot()
		get_path(self, 1, 1)
		snapshot2 = tracemalloc.take_snapshot()
		self.top_stats = snapshot2.compare_to(snapshot1, 'lineno')
		

	def print_results(self, t1, t0):
	
		print("\n\n==========   Needleman–Wunsch  ==========")
		print("Total sequence length is: ", len(self.seq1) + len(self.seq2))
		print(self.path_matrix)
		print('Max {} alignment Score is: {}'.format(self.mode, self.max_score))
		print("Total run time in seconds: ", str(round(t1-t0, 4)))
		count = 1
		if len(self.alignments[0][0]) == 0:
			return
		for alignment in self.alignments:
			print("Alignmnet #", count)
			print(alignment[0])
			print(alignment[1])
			print()
			count += 1
		print("Memory Usage:")
		for stat in self.top_stats[:5]:
			print(stat)
		print()

	def get_all_alignments(self):
	
		return self.alignments
	 

if __name__ == "__main__":
	file = sys.argv[1]
	if not file.endswith('.csv'):
		raise ValueError("Only .csv is accepted for score_matrix") 
	t0 = time.time()
	sequence = SequenceInit("AGTACGGTACGTAA", "TAGAAGTT", pd.read_csv(file), "local")
	t1 = time.time()
	sequence.print_results(t1, t0)