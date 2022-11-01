import sys
import numpy as np


def InputParser(file_path: str) -> dict:
    params = dict()

    with open(file_path, 'r') as file:
        k, _ = list(map(int, file.readline().split()))  # we only need k
        params['k'] = k

    params['data'] = np.loadtxt(file_path, delimiter=' ', skiprows=1)

    return params


def LloydClustering(data: np.array, k: int) -> None:
    centers = np.array([data[i] for i in range(k)])  # initialize random centers

    clusters = [np.inf for _ in range(len(data))]

    while True:
        for i in range(len(data)):
            distances = [np.sqrt(((data[i] - center)**2).sum()) for center in centers]
            clusters[i] = distances.index(min(distances))

        new_centers = np.array([np.mean([data[i] for i in range(len(data)) if clusters[i] == j], axis=0) for j in range(k)])

        if np.all(centers == new_centers):
            break

        centers = new_centers

    for center in centers:
        center = list(map(lambda x: format(round(x, 3), '.3f'), center))
        print(' '.join(list(map(str, center))))


if __name__ == '__main__':
    params = InputParser(sys.argv[1])

    LloydClustering(params['data'], params['k'])