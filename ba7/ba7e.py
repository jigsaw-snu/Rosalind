import sys
import numpy as np


def InputParser(file_path: str) -> dict:
    params = dict()

    with open(file_path, 'r') as file:
        params['n'] = int(file.readline().rstrip())

        params['dist_mat'] = np.zeros(
            shape=(params['n'], params['n']),
            dtype=float
        )

        idx = 0
        while True:
            line = file.readline().rstrip()

            if not line:
                break

            params['dist_mat'][idx] = \
                list(map(int, line.split()))
            
            idx += 1
    
    return params


def MakeTotalDistance(dist_mat: np.array) -> np.array:
    n = dist_mat.shape[0]
    ret = np.zeros(shape=dist_mat.shape)

    # either axis=0 (col-wise) or axis=1 (row-wise) will work fine
    total_dist = dist_mat.sum(axis=1)
    print("<< Total Distance Vector >>", total_dist, '\n', sep='\n')
    for i in range(dist_mat.shape[0]):
        for j in range(dist_mat.shape[0]):
            if i != j:
                ret[i, j] = (n-2)*dist_mat[i, j] - total_dist[i] - total_dist[j]

    print('<< Total Distance >>', ret, '\n', sep='\n')    
    return ret


def GetMinimumDistance(dist_mat: np.array) -> tuple:
    mask = np.zeros(shape=dist_mat.shape, dtype=bool)
    
    for i in range(dist_mat.shape[0]):
        mask[i, i] = True

    masked = np.ma.array(dist_mat, mask=mask)

    print("<< min Index >>", np.unravel_index(np.argmin(masked), shape=dist_mat.shape), '\n', sep='\n')

    return np.unravel_index(np.argmin(masked), shape=dist_mat.shape)


def ReduceLimb(gl_n: int, idx: tuple, dist_mat: np.array) -> list:
    n = dist_mat.shape[0]

    # save limbs
    delta = (dist_mat.sum(axis=1)[idx[0]] - dist_mat.sum(axis=1)[idx[1]]) / (gl_n - 2)
    limbs = np.array([dist_mat[idx[0], idx[j]] + delta, dist_mat[idx[0], idx[1]] - delta]) / 2

    ret = np.zeros(shape=(n, n), dtype=float)

    non_mask = list(range(n))
    non_mask.remove(idx[0])
    non_mask.remove(idx[1])

    for i in range(n):
        for j in range(i+1, n):
            if i in idx or j in idx:
                continue

            ret[i, j] = ret[j, i] = dist_mat[i, j]

    new_node = np.zeros(n-1, dtype=float)

    for i in range(len(non_mask)):
        new_node[i] = (dist_mat[non_mask[i], idx[0]] + dist_mat[non_mask[i], idx[1]] - dist_mat[idx[0], idx[1]]) / 2


    ret = np.delete(ret, idx[1], 1)  # remove high index first to avoid idexing problems
    ret = np.delete(ret, idx[0], 1)
    ret = np.delete(ret, idx[1], 0)
    ret = np.delete(ret, idx[0], 0)

    ret = np.c_[ret, [0]*(n-2)]
    ret = np.r_[ret, [new_node]]  # np.r_[df, [[new_row]]]
    ret = np.delete(ret, n-2, 1)
    ret = np.c_[ret, new_node]  # np.c_[df, [new_col]]
    print("<< Reduced Matrix >>", ret, '\n', sep='\n')
    return ret


def NeighborJoiningTree(gl_n: int, n: int, dist_mat: np.array, red_cnt: int= 0, limbs: dict = dict()) -> list:
    print('\n\n')
    print("<< Dist Matrix >>", dist_mat, '\n', sep='\n')
    print("<< reduce count >>", red_cnt, '\n', sep='\n')
    ret = []

    # Base case
    if n == 2:
        return str(gl_n+red_cnt-2) + '->' + str(gl_n+red_cnt-1) + ':' + str(dist_mat[0, 1])

    new_dist = ReduceLimb(GetMinimumDistance(MakeTotalDistance(dist_mat)), dist_mat)
    red_cnt += 1

    ret.append(NeighborJoiningTree(gl_n, n-1, new_dist, red_cnt))

    return [ret, limbs]


if __name__ == '__main__':
    params = InputParser(sys.argv[1])
    #print(params['n'], params['dist_mat'], sep='\n')
    
    #print(MakeTotalDistance(params['dist_mat']))
    #print(GetMinimumDistance(params['dist_mat']))
    #print(ReduceLimb((0, 1), params['dist_mat']))

    print(NeighborJoiningTree(params['n'], params['n'], params['dist_mat']), sep='\n')