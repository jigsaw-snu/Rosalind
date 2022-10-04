import sys
import numpy as np


def InputParser(file_path: str) -> list:
    ret = []

    with open(file_path, 'r') as file:
        ret += [file.readline().rstrip() for _ in range(2)]

    return ret


def GetEditDistance(seq1: str, seq2: str) -> int:
    score_matrix = np.zeros(shape=(len(seq1)+1, len(seq2)+1))

    for i in range(1, len(seq1)+1):
        for j in range(1, len(seq2)+1):
            mismatch_score = 1
            if seq1[i-1] == seq2[j-1]:
                mismatch_score = 0

            score_matrix[i, j] = min([
                                        score_matrix[i-1, j-1] + mismatch_score,
                                        score_matrix[i-1, j] + 1,
                                        score_matrix[i, j-1] + 1
                                     ])

    return int(score_matrix[len(seq1), len(seq2)])


if __name__ == '__main__':
    seqs = InputParser(sys.argv[1])
    print(GetEditDistance(seqs[0], seqs[1]))
