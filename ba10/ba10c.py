import sys, itertools
import numpy as np


def InputParser(file_path: str) -> dict:
    params = dict()
    
    with open(file_path, 'r') as file:
        params['emission'] = file.readline().rstrip()
        _ = file.readline()  # skip '----'
        
        params['emit_dict'] = dict((y, x) for x, y in enumerate(file.readline().rstrip().split()))
        _ = file.readline()  # skip '----'
        
        params['status_dict'] = dict((y, x) for x, y in enumerate(file.readline().rstrip().split()))
        _ = file.readline()  # skip '----'
        
        _ = file.readline()  # skip header
        params['transition_mat'] = []
        for _ in range(len(params['status_dict'])):  # row is status
            params['transition_mat'].append(list(map(float, file.readline().rstrip().split()[1:])))
        params['transition_mat'] = np.array(params['transition_mat'])
        _ = file.readline()  # skip '----'
        
        _ = file.readline()  # skip header
        params['emit_mat'] = []
        for _ in range(len(params['status_dict'])):  # row is status
            params['emit_mat'].append(list(map(float, file.readline().rstrip().split()[1:])))
        params['emit_mat'] = np.array(params['emit_mat'])
        
    return params


def ViterbiAlgorithm(emission: str, emit2int: dict, status2int: dict, 
                     transition_mat: np.array, emit_mat: np.array) -> str:
    
    score_matrix = np.zeros(shape=(len(status2int), len(emission)))  # for scoring
    direction_matrix = np.zeros(shape=(len(status2int), len(emission)-1), dtype=int)  # for backtracing
    direction_map = dict((y, x) for x, y in enumerate(itertools.product(range(len(status2int)), range(len(status2int)))))  # helper dictionary for direction mark
    #print('Direction map', direction_map, sep='\n')
    #print('\n')
    
    # initialize
    score_matrix[:, 0] = np.log(0.5) + np.log(emit_mat[:, emit2int[emission[0]]])
    
    # forward
    for i in range(1, len(emission)):  # col : emission
        for j in range(len(status2int)):  # row : status; j : current status
            #print('i :', i, ',  j :', j, ',  emit :', emission[i])
            #print("prev scores :", score_matrix[:, i-1])
            #print("trainsition mat[cur_state] :", transition_mat[:, j])
            
            prev_status, max_val = max(enumerate(score_matrix[:, i-1] + np.log(transition_mat[:, j])), key=lambda x: x[1])
            
            #print('prev_status :', list(status2int.keys())[prev_status], ',  max_value :', max_val)
            #print()
            
            score_matrix[j, i] = np.log(emit_mat[j, emit2int[emission[i]]]) + max_val
            direction_matrix[j, i-1] = direction_map[(prev_status, j)]
        
        #print('Column Work Done :', score_matrix[:, i])
        #print('max index :', direction_matrix[:, i-1])
        #print('\n')
    
    #print(score_matrix)
    #print(direction_matrix)
    
    # backtrace
    pointer = len(emission)-2
    end_status_idx = np.argmax(score_matrix[:, pointer])
    hidden_state = list(status2int.keys())[end_status_idx]
    cur_status_idx, _ = list(direction_map.keys())[direction_matrix[end_status_idx, pointer]]
    
    while pointer >= 0:
        pointer -= 1
        hidden_state += list(status2int.keys())[cur_status_idx]
        cur_status_idx, _ = list(direction_map.keys())[direction_matrix[cur_status_idx, pointer]]
    
    
    return hidden_state[::-1]


if __name__ == '__main__':
    params = InputParser(sys.argv[1])
    print(ViterbiAlgorithm(params['emission'], params['emit_dict'], params['status_dict'],
                           params['transition_mat'], params['emit_mat']))