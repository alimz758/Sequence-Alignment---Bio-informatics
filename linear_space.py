# Followed https://en.wikipedia.org/wiki/Hirschberg's_algorithm for this one
import numpy as np
import pandas as pd
import time
import sys
import tracemalloc


COL_MAP = {
	'A': 0,
	'C': 1,
	'T': 2,
	'G': 3,
	'-': 4
}

# get the index of sum max of sum elements
def get_y_mid(scoreL, scoreR):
	max_index = 0
	max_score = float('-Inf')
	for i, (l, r) in enumerate(zip(scoreL, scoreR[::-1])):
		if sum([l, r]) > max_score:
			max_score = sum([l, r])
			max_index = i
	return max_index 

# get the optimal local point; using two rows of the matrix
def get_opt_points(self, seq1, seq2):
	matrix = [[0 for i in range(len(seq2)+1)], [0]]
	max_score = 0
	opt_i = 0
	opt_j = 0
	
	for i in range(1, len(seq1) + 1):
		for j in range(1, len(seq2)+1):
			matrix[1].append(max(matrix[0][j-1] + self.score_matrix.iloc[COL_MAP[seq1[i-1]]][seq2[j-1]], 
                        	matrix[1][j-1] + self.score_matrix.iloc[0]['-'], 
                         	matrix[0][j] + self.score_matrix.iloc[4]['A'],
                          	0))
			if(matrix[1][j] > max_score):
				max_score = matrix[1][j]
				opt_i = i-1
				opt_j = j-1

		matrix[0] = matrix[1][:]
		matrix[1] = [0]
  
	self.max_score = max_score		
	return opt_i, opt_j


# trim sequences for local alignment
def getCropped(self, seq1, seq2):
	opt_i, opt_j = get_opt_points(self, seq1, seq2)
	new_seq1 = seq1[:opt_i+1]
	new_seq2 = seq2[:opt_j+1]
	opt_i, opt_j = get_opt_points(self, new_seq1[::-1], new_seq2[::-1])

	return new_seq1[::-1][:opt_i+1][::-1], new_seq2[::-1][:opt_j+1][::-1]


# gets the last line of the Needleman-Wunsch matrix
def NWScore(self, seq1, seq2):

	len1 = len(seq1) + 1
	len2 = len(seq2) + 1
	last_line = [0] * (len2)
	current_line = [0] * (len2)

	for j in range(1, len2):
		last_line[j] = last_line[j - 1] + self.score_matrix.iloc[4]['A']

	for i in range(1, len1):
		current_line[0] = self.score_matrix.iloc[0]['-'] + last_line[0]
		for j in range(1, len2):
			current_line[j] = max(last_line[j - 1] + self.score_matrix.iloc[COL_MAP[seq1[i-1]]][seq2[j-1]], 
							last_line[j] + self.score_matrix.iloc[0]['-'], 
							current_line[j - 1] + self.score_matrix.iloc[4]['A'])		

		last_line = current_line
		current_line = [0] * (len2)

	return last_line
	
	
# linear space global alignment
def global_alignment(self, seq1, seq2):
	len1 = len(seq1)+1
	len2 = len(seq2)+1
	matrix = np.zeros([len1, len2], dtype='i') 
	
	for i in range(1, len1):
		matrix[i][0] = i * self.score_matrix.iloc[4]['A']
      
	for i in range(1, len2):
		matrix[0][i] = i * self.score_matrix.iloc[0]['-']

	for i in range(1, len1):
		for j in range(1, len2):
			matrix[i][j] = max(matrix[i-1][j-1] + self.score_matrix.iloc[COL_MAP[seq1[i-1]]][seq2[j-1]],
					   matrix[i-1][j] + self.score_matrix.iloc[4]['A'], 
					   matrix[i][j-1] + self.score_matrix.iloc[0]['-'])

	aligned_seq1 = []
	aligned_seq2 = []
	i = len(seq1)
	j = len(seq2)
	while i>0 or j>0:
		if i>0 and  j>0 and matrix[i][j] == matrix[i-1][j-1] + self.score_matrix.iloc[COL_MAP[seq1[i-1]]][seq2[j-1]]:
			aligned_seq1.insert(0, seq1[i-1])			
			aligned_seq2.insert(0, seq2[j-1])
			i -= 1
			j -= 1
		elif i>0 and matrix[i][j] == matrix[i-1][j] + self.score_matrix.iloc[4]['A']:
			aligned_seq1.insert(0, seq1[i-1])
			aligned_seq2.insert(0, '-')
			i -= 1
		elif j>0 and matrix[i][j] == matrix[i][j-1] + self.score_matrix.iloc[0]['-']:
    			aligned_seq1.insert(0, '-')
			aligned_seq2.insert(0, seq2[j-1])
			j -= 1
	aligned_seq1, aligned_seq2 = map(lambda x: "".join(x), [aligned_seq1, aligned_seq2]) 

	return aligned_seq1, aligned_seq2  


def Hirschberg(self, seq1, seq2):
	aligned_seq1 = ""
	aligned_seq2 = ""
	len1 = len(seq1)
	len2 = len(seq2)
 
	if len(seq1) == 0:
		aligned_seq2 = '-' * len2
		aligned_seq1 = seq2
	elif len(seq2) == 0:
		aligned_seq2 = seq1
		aligned_seq1 = '-' * len1
	elif len1 == 1 or len2 == 1:
		aligned_seq1, aligned_seq2 = global_alignment(self, seq1, seq2)
	else:

		xmid = len1 // 2
		scoreL = NWScore(self, seq1[:xmid], seq2)
		scoreR = NWScore(self, seq1[xmid:][::-1], seq2[::-1])
		ymid = get_y_mid(scoreL, scoreR)
		rowLeft, columnUp = Hirschberg(self, seq1[:xmid], seq2[:ymid])
		rowRight, columnDown = Hirschberg(self, seq1[xmid:], seq2[ymid:])
		aligned_seq1 = rowLeft + rowRight
		aligned_seq2 = columnUp + columnDown

	return aligned_seq1, aligned_seq2


class LinearSpaceAlignment():

	def __init__(self, seq1, seq2, score_matrix, mode):
		self.seq1 = seq1
		self.seq2 = seq2
		self.mode = mode
		self.score_matrix = 0
		self.row_size = len(seq1) + 1
		self.col_size = len(seq2) + 1
		self.score_matrix = score_matrix
		tracemalloc.start()
		self.t0 = time.time()
		snapshot1 = tracemalloc.take_snapshot()
		if self.mode == "local":
			self.seq1, self.seq2 = getCropped(self, self.seq1, self.seq2)
		self.aligned_seq1, self.aligned_seq2 = Hirschberg(self, self.seq1, self.seq2)
		self.t1 = time.time()
		snapshot2 = tracemalloc.take_snapshot()
		self.top_stats = snapshot2.compare_to(snapshot1, 'lineno')


	def print_results(self):
		
		print("\n\n==========   Hirschberg  ==========")
		#print('Max {} alignment Score is: {}'.format(self.mode, self.max_score))
		print("Total run time in seconds: ", str(round(self.t1 - self.t0, 4)))
		if len(self.aligned_seq1) == 0:
			return
		print("Alignmnet #", 1)
		print(self.aligned_seq1)
		print(self.aligned_seq2)
		print()
		print("Memory Usage(using `tracemalloc`, the first stats):")
		for stat in self.top_stats[:1]:
			print(str(stat).rsplit(":", 1)[1])
		print()


#test
if __name__ == "__main__":
	file = sys.argv[1]
	if not file.endswith('.csv'):
		raise ValueError("Only .csv is accepted for score_matrix")
	LinearSpaceAlignment("AGTACGGTACGTAA", "TAGAAGTT", pd.read_csv(file), "local")
