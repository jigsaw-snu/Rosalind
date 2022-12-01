import sys
import numpy as np


def InputParser(file_path: str) -> dict:
    params = dict()
    
    with open(file_path, 'r') as file:
        params['emission'] = file.readline().rstrip()
        _ = file.readline()  # skip '----'
        
        params['emit_dict'] = dict((y, x) for x, y in enumerate(file.readline().rstrip().split()))
        _ = file.readline()  # skip '----'
        
        params['path'] = file.readline().rstrip()
        _ = file.readline()  # skip '----'
        
        params['status_dict'] = dict((y, x) for x, y in enumerate(file.readline().rstrip().split()))
        _ = file.readline()  # skip '----'
        
        _ = file.readline()  # skip header
        params['emit_per_status'] = []
        for _ in range(len(params['status_dict'])):
            params['emit_per_status'].append(list(map(float, file.readline().rstrip().split()[1:])))
        params['emit_per_status'] = np.array(params['emit_per_status'])
        
    return params


def GetPrXgivenPi(emission: str, hidden_path: str, 
                  emit2int: dict, status2int: dict, 
                  emit_per_status: np.array) -> float:
    prob = 1.0
    for status, emit in zip(hidden_path, emission):
        prob *= emit_per_status[status2int[status], emit2int[emit]]
    
    return prob


if __name__ == '__main__':
    params = InputParser(sys.argv[1])
    print(GetPrXgivenPi(params['emission'], params['path'], 
                        params['emit_dict'], params['status_dict'],
                        params['emit_per_status']))