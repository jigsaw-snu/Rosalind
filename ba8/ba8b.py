import sys
import numpy as np


def InputParser(file_path: str) -> dict:
    params = dict()
    params['centers'] = list()
    params['data'] = list()

    data_flag = False

    with open(file_path, 'r') as file:
        _ = file.readline()  # we do not need k and m

        while True:
            line = file.readline().rstrip()

            if not line:
                break

            if '-' in line:
                data_flag = True
                continue

            if data_flag:
                params['data'].append(list(map(float, line.split())))
            else:
                params['centers'].append(list(map(float, line.split())))
    
    params['centers'] = np.array(params['centers'])
    params['data'] = np.array(params['data'])

    return params


def SquaredErrorDistortion(centers: np.array, data: np.array) -> float:
    distances = np.array([np.min([np.sqrt(((point - center)**2).sum()) for center in centers]) for point in data])
    distortion = np.mean(distances**2)

    return round(distortion.item(), 3)


if __name__ == '__main__':
    params = InputParser(sys.argv[1])
    
    print(SquaredErrorDistortion(params['centers'], params['data']))