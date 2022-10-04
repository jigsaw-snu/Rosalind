import sys, copy
import numpy as np


INDEL_PENALTY = 3


def InputParser(file_path: str) -> list:
    ret = []

    with open(file_path, 'r') as file:
        ret += [file.readline().rstrip() for _ in range(2)]

    return ret


def GetEditDistance(seq1: str, seq2: str) -> int:
    score_matrix = np.zeros(shape=(len(seq1)+1, len(seq2)+1))
    backtrace = copy.deepcopy(score_matrix)

    for i in range(1, len(seq1)+1):
        for j in range(1, len(seq2)+1):
            MATCH_SCORE = 0
            if seq1[i-1] == seq2[j-1]:
                MATCH_SCORE = 5

            score_matrix[i, j] = max([
                                        score_matrix[i-1, j-1] + MATCH_SCORE,
                                        score_matrix[i-1, j] - INDEL_PENALTY,
                                        score_matrix[i, j-1] - INDEL_PENALTY
                                     ])

            if score_matrix[i, j] == score_matrix[i-1, j-1] + MATCH_SCORE:
                backtrace[i, j] = 3

            elif score_matrix[i, j] == score_matrix[i-1, j] - INDEL_PENALTY:
                backtrace[i, j] = 2

            elif score_matrix[i, j] == score_matrix[i, j-1] - INDEL_PENALTY:
                backtrace[i, j] = 1

    i = len(seq1)
    j = len(seq2)

    ret = ['', '']
    edit_dist = 0

    while i > 0 and j > 0:
        if backtrace[i, j] == 3:
            if seq1[i-1] != seq2[j-1]:
                edit_dist += 1

            ret[0] += seq1[i-1]
            ret[1] += seq2[j-1]
            i -= 1
            j -= 1

        elif backtrace[i, j] == 2:
            edit_dist += 1

            ret[0] += seq1[i-1]
            ret[1] += '-'
            i -= 1

        elif backtrace[i, j] == 1:
            edit_dist += 1

            ret[0] += '-'
            ret[1] += seq2[j-1]
            j -= 1

    while i != 0 or j != 0:
        edit_dist += 1

        if i != 0:
            ret[0] += seq1[i-1]
            i -= 1

        elif j != 0:
            ret[1] += seq2[j-1]
            j -= 1

    return edit_dist


if __name__ == '__main__':
    seqs = InputParser(sys.argv[1])
    print(GetEditDistance(seqs[0], seqs[1]))
