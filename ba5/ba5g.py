import sys, re
import numpy as np


def BlosumParser(blosum_path: str) -> list:
    blosum_ret = [[]]
    blosum_mat = []

    with open(blosum_path, 'r') as blosum:
        while True:
            line = blosum.readline()

            if not line:
                break

            if '#' in line:
                continue

            if not re.findall(r'\d+', line):
                blosum_ret[0] += line.split()
                continue

            blosum_mat.append(list(map(int, line.split()[1:])))

    blosum_ret.append(np.array(blosum_mat))

    # print(blosum_ret[0], blosum_ret[1], sep='\n')

    return blosum_ret


def InputParser(input_path: str) -> list:
    with open(input_path, 'r') as file:
        input_ret = [file.readline().rstrip() for _ in range(2)]

    # print(input_ret[0], input_ret[1], sep='\n')

    return input_ret


def GetEditDistance(seq1: str, seq2: str, sub_mat_info: list, penalty: int) -> int:
    sub_mat_idx = dict((v, k) for (k, v) in enumerate(sub_mat_info[0]))
    sub_mat = sub_mat_info[1]

    score_matrix = np.zeros(shape=(len(seq1)+1, len(seq2)+1))
    score_matrix[:, 0] = [0 - i*penalty for i in range(len(seq1)+1)]
    score_matrix[0, :] = [0 - i*penalty for i in range(len(seq2)+1)]

    backtrace = np.zeros(shape=(len(seq1)+1, len(seq2)+1))

    for i in range(1, len(seq1)+1):
        for j in range(1, len(seq2)+1):
            diag = score_matrix[i-1, j-1] + \
                   sub_mat[sub_mat_idx[seq1[i-1]], sub_mat_idx[seq2[j-1]]]

            down = score_matrix[i-1, j] - penalty

            right = score_matrix[i, j-1] - penalty

            max_score = max([diag, down, right])
            score_matrix[i, j] = max_score

            if max_score == diag:
                backtrace[i, j] = 3

            elif max_score == down:
                backtrace[i, j] = 2

            elif max_score == right:
                backtrace[i, j] = 1

    edit_dist = 0
    aligned = ['', '']

    i = len(seq1)
    j = len(seq2)

    while i > 0 and j > 0:
        if backtrace[i, j] == 3:  # match or mismatch
            if seq1[i-1] != seq2[j-1]:
                edit_dist += 1
            aligned[0] += seq1[i-1]
            aligned[1] += seq2[j-1]
            i -= 1
            j -= 1

        elif backtrace[i, j] == 2:  # down ; gap in seq2
            edit_dist += 1
            aligned[0] += seq1[i-1]
            aligned[1] += '-'
            i -= 1

        elif backtrace[i, j] == 1:  # right : gap in seq1
            edit_dist += 1
            aligned[0] += '-'
            aligned[1] += seq2[j-1]
            j -= 1

    while i != 0 or j != 0:
        if i == 0 and j != 0:
            edit_dist += 1
            aligned[0] += '-'
            aligned[1] += seq2[j-1]
            j -= 1

        elif i != 0 and j == 0:
            edit_dist += 1
            aligned[0] += seq1[i-1]
            aligned[1] += '-'
            i -= 1

    print(aligned[0][::-1], aligned[1][::-1], sep='\n')

    return edit_dist


if __name__ == '__main__':
    input_params = InputParser(sys.argv[1])
    blosum_params = BlosumParser(sys.argv[2])

    print(GetEditDistance(input_params[0], input_params[1], blosum_params, 5))
