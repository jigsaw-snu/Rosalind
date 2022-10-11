import sys
import numpy as np
from itertools import combinations


def InputParser(file_path: str) -> list:
    params = [None, None, None]
    
    with open(file_path, 'r') as file:
        # first line : number of leaf nodes
        params[0] = int(file.readline().rstrip())
        
        # second line : target node (0-based)
        params[1] = int(file.readline().rstrip())
        
        # from third line : Distance Matrix
        params[2] = []
        
        while True:
            line = file.readline().rstrip()
            
            if not line:
                break
            
            params[2].append([int(x) for x in line.split()])
        
        params[2] = np.array(params[2])
        
    return params


def CalculateLimbLength(tree_info: list) -> int:
    pool = list(range(tree_info[0]))
    pool.remove(tree_info[1])
    
    combi = list(combinations(pool, 2))
    
    vals = [(tree_info[2][tree_info[1]][i] + tree_info[2][tree_info[1]][j] - tree_info[2][i][j]) / 2 for (i, j) in combi]
    
    return int(min(vals))


if __name__ == '__main__':
    params = InputParser(sys.argv[1])
    
    #print(params[0], params[1], params[2], sep='\n')
    
    print(CalculateLimbLength(params))