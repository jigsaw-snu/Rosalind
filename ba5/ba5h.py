import sys, copy
import numpy as np


def InputParser(file_path: str) -> list:
    ret = []

    with open(file_path, 'r') as file:
        ret += [file.readline().rstrip() for _ in range(2)]

    return ret


def FittingAlignment(seq1: str, seq2: str) -> list:
    score_matrix = np.zeros(shape=(len(seq1)+1, len(seq2)+1))
    backtrace = copy.deepcopy(score_matrix)

    for i in range(1, len(seq1)+1):
        for j in range(1, len(seq2)+1):
            if seq1[i-1] == seq2[j-1]:
                match_score = 1
            else:
                match_score = -1

            diag = score_matrix[i-1, j-1] + match_score
            down = score_matrix[i-1, j] - 1
            right = score_matrix[i, j-1] - 1

            score_matrix[i, j] = max([0, diag, down, right])

            if score_matrix[i, j] == diag:
                backtrace[i, j] = 3

            elif score_matrix[i, j] == down:
                backtrace[i, j] = 2

            elif score_matrix[i, j] == right:
                backtrace[i, j] = 1

    print(score_matrix, backtrace, sep='\n\n')

    i = len(seq1)
    j = len(seq2)

    ret = [int(score_matrix[i, j]), '', '']

    while j > 0:
        if backtrace[i, j] == 3:



if __name__ == '__main__':
    seqs = InputParser(sys.argv[1])
    fit_aligned = FittingAlignment(seqs[0], seqs[1])

    # print(fit_aligned[0], fit_aligned[1], fit_aligned[2], sep='\n')
