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


def GetRawIndex(idx: int, cur_data: np.array, raw_data: np.array) -> int:
    """Get index from raw data of which item corresponds to item of current data

    Args:
        idx (int): index of item from current data
        cur_data (np.array): current data
        raw_data (np.array): raw data

    Returns:
        int: index of same item from raw data
    """

    return np.where(np.all(raw_data == cur_data[idx], axis=1))[0].item()


def GetDistance(pt1: np.array, pt2: np.array) -> float:
    """Get distance between given point1 and point2

    Args:
        pt1 (np.array): first point
        pt2 (np.array): second point

    Returns:
        float: distance between point1 and point2
    """

    return np.sqrt(np.array([x**2 for x in pt1-pt2]).sum())


def FindMaxPoint(data: np.array, start: int) -> int:
    """Find farthest point from given point

    Args:
        data (np.array): bag of data points
        start (int): index of data point. calculate distance from this point

    Returns:
        int: index of farthest point
    """

    max_idx = 0
    max_dist = 0
    
    for i in range(len(data)):
        dist = GetDistance(data[start], data[i])

        if dist > max_dist:
            max_dist = dist
            max_idx = i
    
    return max_idx


def FarthestFirstTraversal(data: np.array, k: int, dim: int) -> list:
    """Find center points by iteratively searching farthest points

    Args:
        data (np.array): bag of data points
        k (int): number of clusters

    Returns:
        list: data points which are centers of each cluster
    """

    # annotate initial cluster as 0 for all data points
    start_idx = 0  # index of data_raw
    centers = [tuple(map(str, data[start_idx]))]  # result object
    clusters = [start_idx for _ in range(len(data))]
    data_cur = copy.deepcopy(data)

    for _ in range(k-1):
        print('data shape :', data_cur.shape)
        far_idx = FindMaxPoint(data_cur, start_idx)  # index from current data
        far_val = data_cur[far_idx]
        start_val = data_cur[start_idx]
        
        start_idx_raw = GetRawIndex(start_idx, data_cur, data)
        far_idx_raw = GetRawIndex(far_idx, data_cur, data)  # index from raw data
        
        print('start_idx :', start_idx, '    start_val :', start_val)
        print('far_idx :', far_idx, '    far_val :', far_val, '    far_idx_raw :', far_idx_raw)
        #print()
        
        centers.append(tuple(map(str, data_cur[far_idx])))
    
        # reallocate cluster to each data point
        '''
        for i in range(len(data_cur)):
            if GetRawIndex(i, data_cur, data) == 345 or GetRawIndex(i, data_cur, data) == 210:
                print('from far :', GetDistance(far_val, data_cur[i]), '    fram start :', GetDistance(start_val, data_cur[i]))
            if GetDistance(far_val, data_cur[i]) < \
               GetDistance(start_val, data_cur[i]):
                clusters[GetRawIndex(i, data_cur, data)] = far_idx_raw
        '''

        print(data[345], clusters[346], sep=' -> ')
        print(data[210], clusters[211], sep=' -> ')
        print()
        #data_cur = data[[x for x, y in enumerate(clusters) if y == far_idx_raw]]
        data_cur = data[[x for x in range(len(data)) if x != start_idx_raw]]
        
        # must recalculate index of previous farthest point
        start_idx = GetRawIndex(far_idx_raw, data, data_cur)

    # printing
    for i in range(k):
        print(' '.join(centers[i]))
    
    return centers


if __name__ == '__main__':
    params = InputParser(sys.argv[1])
    #print(params['k'], params['dim'])
    #print(params['data'])
    
    _ = FarthestFirstTraversal(params['data'], params['k'], params['dim'])