import sys, copy
import numpy as np


MATCH_SCORE = 1
INDEL_PENALTY = 1

def InputParser(file_path: str) -> list:
    ret = []

    with open(file_path, 'r') as file:
        ret += [file.readline().rstrip() for _ in range(2)]

    return ret


def FittingAlignment(seq1: str, seq2: str) -> list:
    score_matrix = np.zeros(shape=(len(seq1)+1, len(seq2)+1))
    backtrace = copy.deepcopy(score_matrix)

    # let len(seq1) > len(seq2)
    # initialize column with 0 --> setting first column to  0 means alignment can be started within any place of seq1
    # do not initialize row with 0 --> setting first row to 0 means alignment can be started within any place of seq2
    # since, we are trying to align full sequence of seq2 to subsequence of seq1
    # --> give indel penalty for seq2
    score_matrix[0, :] = [-i*INDEL_PENALTY for i in range(len(seq2)+1)]

    max_score = 0
    max_index = [0, 0]

    for i in range(1, len(seq1)+1):
        for j in range(1, len(seq2)+1):
            if seq1[i-1] != seq2[j-1]:
                MATCH_SCORE = -1

            # use global alignment in fitting alignment
            # we want to align whole sequence of seq2 to seq1
            diag = score_matrix[i-1, j-1] + MATCH_SCORE
            down = score_matrix[i-1, j] - INDEL_PENALTY
            right = score_matrix[i, j-1] - INDEL_PENALTY

            score_matrix[i, j] = max([diag, down, right])

            if j == len(seq2):
                if max_score <= score_matrix[i, j]:
                    max_score = score_matrix[i, j]
                    max_index = [i, j]

            if score_matrix[i, j] == diag:
                backtrace[i, j] = 3

            elif score_matrix[i, j] == down:
                backtrace[i, j] = 2

            elif score_matrix[i, j] == right:
                backtrace[i, j] = 1

    print(score_matrix, backtrace, '\n', sep='\n\n')

    # start backtracking from index of maximum value in last column
    # starting backtracing from last column means starting backtracking from last sequence of seq2
    i, j = max_index
    ret = [int(max_score), '', '']

    while i > 0 and j > 0:
        if backtrace[i, j] == 3:
            ret[1] += seq1[i-1]
            ret[2] += seq2[j-1]
            i -= 1
            j -= 1

        elif backtrace[i, j] == 2:
            ret[1] += seq1[i-1]
            ret[2] += '-'
            i -= 1

        elif backtrace[i, j] == 1:
            ret[1] += '-'
            ret[2] += seq2[j-1]
            j -= 1

    # all sequence of seq2 must be aligned! be sure to check j != 0
    if j > 0:
        ret[1] += '-'
        ret[2] += seq2[j-1]

    print(ret[0], ret[1], ret[2], sep='\n')


if __name__ == '__main__':
    seqs = InputParser(sys.argv[1])
    fit_aligned = FittingAlignment(seqs[0], seqs[1])

    # print(fit_aligned[0], fit_aligned[1], fit_aligned[2], sep='\n')
