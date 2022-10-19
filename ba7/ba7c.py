import sys
from itertools import combinations
import numpy as np


def InputParser(file_path: str) -> dict:
    params = dict()

    with open(file_path, 'r') as file:
        params['n'] = int(file.readline().rstrip())
        params['dist_mat'] = np.zeros(shape=(params['n'], params['n']), dtype=int)

        mat_idx = 0
        while True:
            line = file.readline().rstrip()

            if not line:
                break

            params['dist_mat'][mat_idx] = list(map(int, line.split()))
            mat_idx += 1
    
    return params


def LimbLength(node: int, dist_mat: np.array) -> int:
    return int(min(
        [(dist_mat[node, i] + dist_mat[node, j] - dist_mat[i, j]) / 2 \
         for i, j in combinations(range(node), 2) if i != node and j != node]
    ))


def TreeBuilder(n: int, dist_mat: np.array) -> str:
    # Base Case : n == 2
    if n == 2:
        return (0, 1, dist_mat[0, 1])

    limb = LimbLength(n, dist_mat)

    # Make Bald Tree : subtract limb length from n-th row and column
    for i in range(1, n):
        dist_mat[i, n] -= limb
        dist_mat[n, i] = dist_mat[i, n]

    # Get limb length from Bald Tree
    ()

    # Make Trimmed Tree : remove n-th row and column



if __name__ == "__main__":
    params = InputParser(sys.argv[1])
    #print(params['n'], params['dist_mat'], sep='\n')

    #print(LimbLength(2, params['dist_mat']))