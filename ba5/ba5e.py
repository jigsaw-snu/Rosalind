import sys, re
import numpy as np


'''
    <BlosumParser>
    Read Blosum matrix in txt form and save into dataframe form
'''
def BlosumParser(blosum_path: str) -> list:
    contents = []
    idx = []

    with open(blosum_path, 'r') as blosum:
        while True:
            line = blosum.readline().strip()

            if '#' in line:
                continue

            if not line:
                break

            if not re.findall(r'\d+', line):
                idx += line.split()
                continue

            contents.append(list(map(int, line.split()[1:])))

    return [idx, np.array(contents)]


'''
    <InputParser>
    Parse input data and save it to dictionary
'''
def InputParser(file_path: str) -> dict:
    params = dict()
    params['sequence'] = []

    with open(file_path, 'r') as file:
        while True:
            line = file.readline().rstrip()

            if not line:
                break

            params['sequence'].append(line)

    return params


'''
    <GlobalAligner>
    Find global alignment between two given sequences
'''
def GlobalAligner(sequence: list, blosum_info: list, penalty: int) -> list:
    blosum = blosum_info[1]  # BLOSUM matrix
    b_idx = dict((v, k) for (k, v) in enumerate(blosum_info[0]))

    score_matrix = np.zeros(shape=(len(sequence[0])+1, len(sequence[1])+1))

    # Don't' forget to initialize first row and first column  with gap penalty
    score_matrix[0, :] = [0 - i*penalty for i in range(len(sequence[1])+1)]
    score_matrix[:, 0] = [0 - i*penalty for i in range(len(sequence[0])+1)]

    # print(score_matrix, '\n')

    # make sure that you need to store backtracking info in matrix not a list
    backtrack = np.zeros(shape=(len(sequence[0])+1, len(sequence[1])+1))

    for i in range(1, len(sequence[0])+1):
        for j in range(1, len(sequence[1])+1):
            max_score = max([
                score_matrix[i-1, j-1] +
                blosum[b_idx[sequence[0][i-1]], b_idx[sequence[1][j-1]]],

                score_matrix[i-1, j] - penalty,

                score_matrix[i, j-1] - penalty
            ])

            score_matrix[i, j] = max_score

            if max_score == score_matrix[i-1, j-1] + \
                    blosum[b_idx[sequence[0][i-1]], b_idx[sequence[1][j-1]]]:
                backtrack[i, j] = 3  # let 3 as diagonal movement
            elif max_score == score_matrix[i-1, j] - penalty:
                backtrack[i, j] = 2  # let 2 as down movement
            elif max_score == score_matrix[i, j-1] - penalty:
                backtrack[i, j] = 1  # let 1 as right movement

            score_matrix[i, j] = max([
                score_matrix[i-1, j-1] +
                blosum[b_idx[sequence[0][i-1]], b_idx[sequence[1][j-1]]],

                score_matrix[i-1, j] - penalty,

                score_matrix[i, j-1] - penalty
            ])

    # print(score_matrix, '\n')
    # print(backtrack, '\n')

    i = len(sequence[0])
    j = len(sequence[1])
    res = [int(score_matrix[i, j]), '', '']

    # Backtracking
    # remember sequence[0] is row, and sequence[1] is column
    while i >= 0 and j >= 0:
        # Handling corner case
        if not i and j:  # i is 0
            res[1] += '-'
            res[2] += sequence[1][j-1]
            break

        elif i and not j:  # j is 0
            res[1] += sequence[0][i-1]
            res[2] += '-'
            break

        elif not i and not j:  # if both i and j are 0, end this loop
            break

        if backtrack[i, j] == 3:  # diagonal : match or mismatch
            res[1] += sequence[0][i-1]
            res[2] += sequence[1][j-1]
            i -= 1
            j -= 1

        elif backtrack[i, j] == 2:  # down : gap in sequence[1]
            res[1] += sequence[0][i-1]
            res[2] += '-'
            i -= 1

        elif backtrack[i, j] == 1:  # right : gap in sequence[0]
            res[1] += '-'
            res[2] += sequence[1][j-1]
            j -= 1

    res[1] = res[1][::-1]
    res[2] = res[2][::-1]

    return res


if __name__ == '__main__':
    blosum = BlosumParser(sys.argv[2])
    params = InputParser(sys.argv[1])
    align = GlobalAligner(params['sequence'], blosum_info=blosum, penalty=5)

    # print(params['sequence'][0], params['sequence'][1], sep='\n')
    print(align[0], align[1], align[2], sep='\n')
