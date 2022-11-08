import sys
import numpy as np


def InputParser(file_path: str) -> dict():
    params = dict()
    
    with open(file_path, 'r') as file:
        params['k'], params['m'] = list(map(int, file.readline().split()))
        params['beta'] = float(file.readline().rstrip())
        
        params['data'] = []
        while True:
            line = file.readline().rstrip()
            
            if not line:
                break
            
            params['data'].append(list(map(float, line.split())))
    
        params['data'] = np.array(params['data'])
    
    return params


def GetDistance(point: np.array, center: np.array) -> float:
    return np.sqrt(np.square(point - center).sum())


def SoftKmeansClustering(k: int, beta: float, datapoints: np.array) -> None:
    # initialize centers with first k points
    centers = datapoints[:k]

    steps = 100
    while steps > 0:
        steps -= 1

        # E-step (Estimating Hidden Matrix)
        hidden_matrix = np.array([[np.exp(-beta*GetDistance(point, center)) for point in datapoints] for center in centers])
        hidden_matrix = hidden_matrix / hidden_matrix.sum(axis=0)  # axis=0 : column-wise
        
        # M-step (Estimating Parameters)
        new_centers = ((hidden_matrix @ datapoints).T / hidden_matrix.sum(axis=1)).T  # axis=1 : row-wise
        
        if np.array_equal(centers, new_centers):
            break
        
        centers = new_centers
    
    for center in centers:
        print(' '.join(list(map(lambda x: format(x, '.3f'), center))))
        
        

if __name__ == '__main__':
    params = InputParser(sys.argv[1])
    SoftKmeansClustering(params['k'], params['beta'], params['data'])