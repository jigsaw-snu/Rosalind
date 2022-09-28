import sys, re
import pandas as pd
import numpy as np


'''
    <BlosumParser>
    Read Blosum matrix in txt form and save into dataframe form
'''
def BlosumParser(blosum_path: str) -> pd.DataFrame:
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

    res = pd.DataFrame(contents, columns=idx, index=idx)

    return res


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
def GlobalAligner(sequence: list, blosum: pd.DataFrame, penalty: int) -> list:
    score_matrix = np.zeros(shape=(len(sequence[0])+1, len(sequence[1])+1))
    print(score_matrix)

    # first row and first column is zero padded
    for i in range(1, len(sequence[0])+1):
        for j in range(1, len(sequence[1])+1):
            score_matrix[i, j] = max([
                score_matrix[i-1, j-1] + \
                    blosum[sequence[0][i-1]][sequence[1][j-1]],
                score_matrix[i-1, j] - penalty,
                score_matrix[i, j-1] - penalty
            ])

    print(score_matrix)

if __name__ == '__main__':
    blosum = BlosumParser(sys.argv[2])
    params = InputParser(sys.argv[1])
    align = GlobalAligner(params['sequence'], blosum = blosum, penalty=5)

   # print(align[0], align[1], sep='\n')
