import sys
import numpy as np


def InputParser(filepath: str) -> dict:
    params = dict()
    
    with open(filepath, 'r') as file:
        params['path'] = file.readline().rstrip()
        _ = file.readline()  # skip '----'
        
        params['dict'] = dict((y, x) for x, y in enumerate(file.readline().rstrip().split()))
        _ = file.readline()  # skip '----'
        
        params['trans_mat'] = []
        _ = file.readline()  # skip header
        for _ in range(len(params['dict'])):
            params['trans_mat'].append(list(map(float, file.readline().rstrip().split()[1:])))
        params['trans_mat'] = np.array(params['trans_mat'])
    
    return params


def GetPathProbability(path: str, chr2int: dict, trans_mat: np.array) -> float:
    prob = 1.0
    prev = path[0]

    for chr in path[1:]:
        prob *= trans_mat[chr2int[prev], chr2int[chr]]
        prev = chr
    
    return prob * 0.5  # initial transition probability = 0.5


if __name__ == '__main__':
    params = InputParser(sys.argv[1])
    print(GetPathProbability(params['path'], params['dict'], params['trans_mat']))