import sys
import pandas as pd
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

    # make zero padded matrix
    # zero padding is needed to get mat[i-1][j-1] without any errors
    mat = pd.DataFrame(0, index=range(n+1), columns=range(m+1))
    backtrack = pd.DataFrame('', index=range(n), columns=range(m))

    # cautious indexing : df[column][row]
    for i in tqdm(range(1, m+1)):
        for j in range(1, n+1):
            if str1[j-1] == str2[i-1]:
                mat[i][j] = mat[i-1][j-1]+1
                backtrack[i-1][j-1] = 'diag'

            elif mat[i][j-1] > mat[i-1][j]:
                mat[i][j] = mat[i][j-1]
                backtrack[i-1][j-1] = 'down'

            elif mat[i][j-1] <= mat[i-1][j]:
                mat[i][j] = mat[i-1][j]
                backtrack[i-1][j-1] = 'right'

    ''' # for debugging purpose
    mat2 = mat.copy()
    mat2.index = list('0'+str1)
    mat2.columns = list('0'+str2)
    print(mat2, '\n')

    backtrack2 = backtrack.copy()
    backtrack2.index = list(str1)
    backtrack2.columns = list(str2)
    print(backtrack2, '\n')
    '''

    LCS = ''
    i, j = m, n
    while j > 0 and i > 0:
        if mat[i][j] == mat[i-1][j-1]+1 and backtrack[i-1][j-1] == 'diag':
            LCS += str2[i-1]
            # go diagonally
            i -= 1
            j -= 1
            continue

        elif mat[i][j] == mat[i-1][j] and backtrack[i-1][j-1] == 'right':
            i -= 1  # go left
            continue

        elif mat[i][j] == mat[i][j-1] and backtrack[i-1][j-1] == 'down':
            j -= 1  # go up
            continue

    return LCS[::-1]


if __name__ == '__main__':
    params = InputParser(sys.argv[1])
    print(len(params['sequence'][0]), len(params['sequence'][1]))
    print(LCSBuilder(params['sequence'][0], params['sequence'][1]))
