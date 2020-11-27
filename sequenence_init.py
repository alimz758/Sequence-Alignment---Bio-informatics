class SequenceInit(object):
    def __init__(self, seq1, seq2, score_matrix):
        self.seq1 = seq1
        self.seq2 = seq2
        self.row_size = len(seq1) + 1
        self.col_size = len(seq2) + 1
        self.score_matrix = score_matrix
        self.indel = self.score_matrix['-']