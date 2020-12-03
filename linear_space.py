# Followed https://en.wikipedia.org/wiki/Hirschberg's_algorithm for this one
import copy
import pdb
import numpy as np
import pandas as pd
import time
import resource
import sys
import tracemalloc


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



def NWScore(self, seq1, seq2):

	len1 = len(seq1) + 1
	len2 = len(seq2) + 1
	pre_row = [0] * (len2)
	cur_row = [0] * (len2)

	for j in range(1, len2):
		pre_row[j] = pre_row[j - 1] + self.score_matrix.iloc[4]['A']

	for i in range(1, len1):
		cur_row[0] = self.score_matrix.iloc[0]['-'] + pre_row[0]
		for j in range(1, len2):
			cur_row[j] = max(pre_row[j - 1] + self.score_matrix.iloc[COL_MAP[seq1[i-1]]][seq2[j-1]], 
							pre_row[j] + self.score_matrix.iloc[0]['-'], 
							cur_row[j - 1] + self.score_matrix.iloc[4]['A'])

		pre_row = cur_row
		cur_row = [0] * (len2)

	return pre_row

def local_alignment(self, seq1, seq2):
	column = len(seq1)+1
	row = len(seq2)+1
	matrix = np.zeros([row, column], dtype='i,O') 
	max_score = float('-inf')
	opt_i = 0
	opt_j = 0
	
	for i in range(1, row):
		for j in range(1, column):
			diag = matrix[i-1][j-1][0] + self.score_matrix.iloc[COL_MAP[self.seq1[j-1]]][self.seq2[i-1]]
			up = matrix[i-1][j][0] + self.score_matrix.iloc[4]['A']
			left = matrix[i][j-1][0] + self.score_matrix.iloc[0]['-']
			matrix[i][j][0] = max(diag, up, left, 0)

			if matrix[i][j][0]==left:
				matrix[i][j][1] = 'L'

			if matrix[i][j][0]==diag:
				matrix[i][j][1] = 'D'

			if matrix[i][j][0]==up:
				matrix[i][j][1] = 'U'
			#for local
			if matrix[i][j][0] >= max_score:
				max_score = matrix[i][j][0]
				opt_i = i
				opt_j = j
	row = []
	column = []
	i = opt_i
	j = opt_j
	while matrix[i][j][1]:
		if matrix[i][j][0] == 0:
			break
		elif matrix[i][j][1] == 'D':
			column.insert(0, seq2[i-1])
			row.insert(0, seq1[j-1])
			i -= 1
			j -= 1
		elif matrix[i][j][1] == 'U':
			row.insert(0, '-')
			column.insert(0, seq2[i-1])
			i -= 1
		elif matrix[i][j][1] == 'L':
			column.insert(0, '-')
			row.insert(0, seq1[j-1])
			j -= 1
	
	row, column = map(lambda x: "".join(x), [row, column]) 
	return column, row 
	
	
def global_alignment(self, seq1, seq2):
	row = len(seq1)+1
	column = len(seq2)+1
	matrix = np.zeros([row, column], dtype='i,O') 
	
	for i in range(1, column):
		matrix[0][i][0] = i * self.score_matrix.iloc[0]['-']
		matrix[0][i][1] = 'L'

	for i in range(1, row):
		matrix[i][0][0] = i * self.score_matrix.iloc[4]['A']
		matrix[i][0][1] = 'U'

	for i in range(1, row):
		for j in range(1, column):
			matrix[i][j][0] = max(matrix[i-1][j-1][0] + self.score_matrix.iloc[COL_MAP[seq1[i-1]]][seq2[j-1]], matrix[i-1][j][0] + self.score_matrix.iloc[4]['A'], matrix[i][j-1][0] + self.score_matrix.iloc[0]['-'])
			if matrix[i][j][0] == matrix[i-1][j-1][0] + self.score_matrix.iloc[COL_MAP[seq1[i-1]]][seq2[j-1]]:
				matrix[i][j][1] =  'D'
			elif matrix[i][j][0] == matrix[i-1][j][0] + self.score_matrix.iloc[4]['A']:
				matrix[i][j][1] = 'U'
			else:
				matrix[i][j][1] = 'L'

	row = []
	column = []
	i = len(seq1)
	j = len(seq2)
	while matrix[i][j][1]:
		if matrix[i][j][1] == 'D':
			column.insert(0, seq2[j-1])
			row.insert(0, seq1[i-1])
			i -= 1
			j -= 1
		elif matrix[i][j][1] == 'U':
			column.insert(0, '-')
			row.insert(0, seq1[i-1])
			i -= 1
		elif matrix[i][j][1] == 'L':
			row.insert(0, '-')
			column.insert(0, seq2[j-1])
			j -= 1
	row, column = map(lambda x: "".join(x), [row, column]) 

	return row, column 


def Hirschberg(self, seq1, seq2):
	row = ""
	column = ""
	len1 = len(seq1)
	len2 = len(seq2)
 
	if len(seq1) == 0:
		column = '-' * len2
		row = seq2
	elif len(seq2) == 0:
		column = seq1
		row = '-' * len1
	elif len1 == 1 or len2 == 1:
		if self.mode == "global":
			row, column = global_alignment(self, seq1, seq2)
		else:
			row, column = local_alignment(self, seq1, seq2)
	else:

		xmid = len1 // 2
		scoreL = NWScore(self, seq1[:xmid], seq2)
		scoreR = NWScore(self, seq1[xmid:][::-1], seq2[::-1])
		ymid = get_max_score(scoreL, scoreR)

		rowLeft, columnUp = Hirschberg(self, seq1[:xmid], seq2[:ymid])
		rowRight, columnDown = Hirschberg(self, seq1[xmid:], seq2[ymid:])
		row = rowLeft + rowRight
		column = columnUp + columnDown
	return row, column


class LinearSpaceAlignment():

	def __init__(self, seq1, seq2, score_matrix, mode):
		self.seq1 = seq1
		self.seq2 = seq2
		self.mode = mode
		self.row_size = len(seq1) + 1
		self.col_size = len(seq2) + 1
		self.score_matrix = score_matrix
		tracemalloc.start()
		self.t0 = time.time()
		snapshot1 = tracemalloc.take_snapshot()
		self.aligned_seq2, self.aligned_seq1 = Hirschberg(self, self.seq1, self.seq2)
		self.t1 = time.time()
		snapshot2 = tracemalloc.take_snapshot()
		self.top_stats = snapshot2.compare_to(snapshot1, 'lineno')


	def print_results(self):
    	
		print("\n\n==========   Hirschberg  ==========")
		print("Total sequence length is: ", len(self.seq1) + len(self.seq2))
		print("Total run time in seconds: ", str(round(self.t1 - self.t0, 4)))
		if len(self.aligned_seq1) == 0:
			return
		print("Alignmnet #", 1)
		print(self.aligned_seq1)
		print(self.aligned_seq2)
		print()
		print("Memory Usage:")
		for stat in self.top_stats[:5]:
			print(stat)
		print()


#test
if __name__ == "__main__":
	file = sys.argv[1]
	if not file.endswith('.csv'):
		raise ValueError("Only .csv is accepted for score_matrix")
	LinearSpaceAlignment("AGTACGGTACGTAA", "TAGAAGTT", pd.read_csv(file), "local")
