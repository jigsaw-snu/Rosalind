import numpy as np


'''
    <Class : Aligner>
    - Substitution Matrix (given)
    - Score Matrix
    - Backtrace Matrix
'''
class Aligner:
    def __init__(self, seq1: str, seq2: str, sub_mat: np.array):
        self.seq1, self.seq2 = seq1, seq2
        self.sub_mat = sub_mat


    def GlobalAligner():
        pass


    def LocalAligner():
        pass

