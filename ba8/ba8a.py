import sys
import numpy as np


def InputParser(file_path: str) -> dict:
    params = dict()

    with open(file_path, 'r') as file:
        params['k'], _ = list(map(int, file.readline().split()))  # we only need k
    
    params['data'] = np.loadtxt(file_path, delimiter=' ', skiprows=1)

    return params


def FarthestFirstTravversal(data: np.array, k: int) -> None:
    centers = [data[0]]  # initialize with first point

    while len(centers) < k:
        distances = [np.min([np.sqrt(((point - center)**2).sum()) for center in centers]) for point in data]
        max_idx = distances.index(max(distances))
        centers.append(data[max_idx])
    
    for center in centers:
        print(' '.join(list(map(str, center))))


if __name__ == '__main__':
    params = InputParser(sys.argv[1])
    #print(params['k'])
    #print(params['data'])

    FarthestFirstTravversal(params['data'], params['k'])