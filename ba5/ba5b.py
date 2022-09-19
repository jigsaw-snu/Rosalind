import sys
import pandas as pd


'''
    <InputParser>
    Parse input data and save it into dictionary
'''
def InputParser(file_path: str) -> dict:
    params = dict()

    with open(file_path, 'r') as file:
        # first line implies n and m
        params['n'], params['m'] = list(map(int, file.readline().rstrip().split()))

        # range(2, '-') : n x (m+1) matrix for down
        # range('-'+1, -1) : (n+1) x m matrix for right
        mat = []

        while True:
            line = file.readline().rstrip()

            if not line:
                params['right'] = pd.DataFrame(mat)
                break

            # '-' separates down weight matrix and right weight matrix
            if '-' in line:
                params['down'] = pd.DataFrame(mat)
                mat = []
                continue

            mat.append(list(map(int, line.split())))

    return params


'''
    <ManhattanTourist>
    Calculate longest path from given wight matrices
'''
def ManhattanTourist(n: int, m: int, Down: pd.DataFrame, Right: pd.DataFrame) -> int:
    mat = pd.DataFrame(0, index=range(n+1), columns=range(m+1))

    for i in range(1, m+1):
        mat[i][0] = mat[i-1][0] + Right[i-1][0]

    for j in range(1, n+1):
        mat[0][j] = mat[0][j-1] + Down[0][j-1]

    for i in range(1, m+1):
        for j in range(1, n+1):
            mat[i][j] = max([mat[i-1][j] + Right[i-1][j], mat[i][j-1] + Down[i][j-1]])

    return mat[m][n]


if __name__ == '__main__':
    params = InputParser(sys.argv[1])

    '''# InputParser test
    print(params['n'], params['m'])
    print(params['down'])
    print(params['right'])
    '''

    print(ManhattanTourist(params['n'], params['m'], params['down'], params['right']))
