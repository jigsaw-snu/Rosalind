import sys
import numpy as np
from tqdm import tqdm


'''
    <InputParser>
    Parse input file and save data into dictionary
'''
def InputParser(file_path: str) -> dict:
    params = dict()

    with open(file_path, 'r') as file:
        params['sequence'] = list(map(str.rstrip, file.readlines()))

    return params


'''
    <LCSBuilder>
    Build LCS from given strings and return list of backtracking information
'''
def LCSBuilder(str1: str, str2: str) -> str:
    n = len(str1)  # row ; index j
    m = len(str2)  # column ; index i

    # make zero padded matrix : shape = (n+1) x (m+1)
    # zero padding is needed to get mat[i-1][j-1] without any errors
    mat = np.zeros((n+1, m+1))
    backtrack = np.zeros((n, m))

    # cautious indexing : df[column][row] in pandas, arr[row, column] in numpy
    for i in tqdm(range(1, m+1)):
        for j in range(1, n+1):
            if str1[j-1] == str2[i-1]:
                mat[j, i] = mat[j-1, i-1]+1
                backtrack[j-1, i-1] = 3  # we encode 3 as diagonal movement

            elif mat[j-1, i] > mat[j, i-1]:
                mat[j, i] = mat[j-1, i]
                backtrack[j-1, i-1] = 2  # we encode 2 as down movement

            elif mat[j-1, i] <= mat[j, i-1]:
                mat[j, i] = mat[j, i-1]
                backtrack[j-1, i-1] = 1  # we encode 1 as right movement

    LCS = ''
    i, j = m, n
    while j > 0 and i > 0:
        '''
        print('i :', i, ', j :', j)
        print('mat[i][j] :', mat[i][j])
        print('backtrack[i-1][j-1] :', backtrack[i-1][j-1])
        print('LCS :', LCS[::-1], '\n')
        '''
        if mat[j, i] == mat[j-1, i-1]+1 and backtrack[j-1, i-1] == 3:
            LCS += str2[i-1]
            # go diagonally
            i -= 1
            j -= 1
            continue

        elif mat[j, i] == mat[j, i-1] and backtrack[j-1, i-1] == 1:
            i -= 1  # go left
            continue

        elif mat[j, i] == mat[j-1, i] and backtrack[j-1, i-1] == 2:
            j -= 1  # go up
            continue

    return LCS[::-1]


if __name__ == '__main__':
    params = InputParser(sys.argv[1])
    print(LCSBuilder(params['sequence'][0], params['sequence'][1]))
