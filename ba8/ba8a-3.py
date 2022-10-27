import sys, copy
import numpy as np

np.set_printoptions(formatter={'float': '{: 0.1f}'.format})


def InputParser(file_path: str) -> dict:
    """Parse input file and return it as a dictionary

    Args:
        file_path (str): path for input data file

    Returns:
        dict: object containing number of cluster k, dimension of data, and data itself
    """
    
    params = dict()
    
    with open(file_path, 'r') as file:
        params['k'], params['dim'] = list(map(int, file.readline().rstrip().split()))

        params['data'] = list()
        
        while True:
            line = list(map(float, file.readline().rstrip().split()))
            
            if not line:
                break
            
            params['data'].append(line)
        
        params['data'] = np.array(params['data'])
    
    return params


def GetDistance(pt1: np.array, pt2: np.array) -> float:
    return np.sqrt(np.array([x**2 for x in pt1-pt2]).sum())


def FarthestFirstTraversal(data: np.array, k: int) -> None:
    centers = [data[0].tolist()]  # initialize centers with first entry
    
    while len(centers) < k:
        # save longest distance from any of centers
        distances = [0 for _ in range(len(data))]
        
        for i in range(len(centers), len(data)):
            if data[i].tolist() in centers:
                continue
            
            for center in centers:
                dist = max([GetDistance(data[i], center)])
                
                if dist > distances[i]:
                    distances[i] = dist
                    
        # data point with maximum distance from any of centers will become next center
        next_cidx = distances.index(max(distances))
        centers.append(data[next_cidx].tolist())
    
    # print results
    for center in centers:
        center = list(map(str, center))
        print(' '.join(center))


if __name__ == '__main__':
    params = InputParser(sys.argv[1])
    #print(params['k'], params['dim'])
    #print(params['data'])

    print('-'*30, 'answer', '-'*30)
    #GetDistance(np.array(0.8, 12.0, 17.5, 0.9, 7.2), np.array(23.1, 31.1, 3.6, 0.8, 0.3))
    #GetDistance(np.array(0.3, 16.4, 8.9, 34.6, 24.6), np.array(23.1, 31.1, 3.6, 0.8, 0.3))
    print(GetDistance(np.array([0.8, 12.0, 17.5, 0.9, 7.2]), np.array([32.3, 1.9, 5.1, 16.2, 8.8])))
    print(GetDistance(np.array([0.3, 16.4, 8.9, 34.6, 24.6]), np.array([32.3, 1.9, 5.1, 16.2, 8.8])))
    print()
    
    print('-'*30, 'my result', '-'*30)
    print(GetDistance(np.array([0.8, 12.0, 17.5, 0.9, 7.2]), np.array([31.1, 2.1, 12.5, 1.1, 2.5])))
    print(GetDistance(np.array([0.3, 16.4, 8.9, 34.6, 24.6]), np.array([31.1, 2.1, 12.5, 1.1, 2.5])))
    #GetDistance(np.array(0.8, 12.0, 17.5, 0.9, 7.2), np.array(32.3, 1.9, 5.1, 16.2, 8.8))
    #GetDistance(np.array(0.3, 16.4, 8.9, 34.6, 24.6), np.array(32.3, 1.9, 5.1, 16.2, 8.8))
    print()
    
    FarthestFirstTraversal(params['data'], params['k'])