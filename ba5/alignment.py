import numpy as np


'''
    <Class : Aligner>
    @variables:
        - seq1, seq2 : target sequences
        - sub_matrix : substitution matrix given as parameter
        - penalty : gap (indel) penalty
        - score_matrix
        - backtrace

    @functions:
        - GlobalAligner
        - LocalAligner
'''
class Aligner:
    def __init__(self, seq1: str, seq2: str, sub_mat: np.array, penalty: int):
        # prepare alignment
        self.seq1, self.seq2 = seq1, seq2
        self.sub_matrix = sub_mat
        self.penalty = penalty
        self.score_matrix = self.backtrace = \
                np.zeros(shape=(len(self.seq1)+1, len(self.seq2)+1))


    def Backtrack(self):
        pass


    def GlobalAligner(self):
        # global alignment must initialize first row/col with gap penalty
        self.score_matrix[0, :] = [0 - i*self.penalty for i in range(len(seq2)+1)]
        self.score_matrix[:, 0] = [0 - i*self.penalty for i in range(len(seq1)+1)]

        # fill out score matrix
        for i in range(1, len(seq1)+1):
            for j in range(1, len(seq2)+1):
                self.score_matrix[i, j] = \
                        


    def LocalAligner(self):
        # fill out score matrix
