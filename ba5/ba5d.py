import re, sys
import numpy as np
from tqdm import tqdm


class Node:
    def __init__(self):
        self.max_val = 0
        self.innode = ''


'''
    <InputParser>
    Parse input data and save it into a dictionary
'''
def InputParser(file_path: str) -> dict:
    params = dict()

    with open(file_path, 'r') as file:
        # first and second line refers to source and sink node
        line = file.readline().rstrip()
        params['source'] = int(line)

        line = file.readline().rstrip()
        params['sink'] = int(line)

        # from third line, it refers to weight information
        # for the ease of access, let's store this data in numpy array
        #
        # Assumptions
        # 1. node[i] is not linked to node[j] for all node[i] > node[j]
        #    --> only smaller node to bigger node is allowed
        # 2. according to 1), let us consider sink node as the maximum node
        params['weight'] = []
        params['minimap'] = np.zeros(shape=(params['sink']+1, params['sink']+1))

        while True:
            line = file.readline().rstrip()

            if not line:
                break

            src, dst, wt = list(map(int, re.findall(r'\d', line)))
            params['weight'].append((src, dst, wt))

            if src <= params['sink'] and dst <= params['sink']:
                params['minimap'][src, dst] = wt
            else:
                continue

    return params


'''
    <PathSelector>
    Select Best Path from source node to sink node
'''
def PathSelector(src: int, dst: int, weight: list) -> list:
    # add different Node objects to nodes list
    # make list index same as iterator below
    nodes = []
    for _ in range(dst+1):
        nodes.append(Node())

    for i in range(src+1, dst+1):
        print('node_vals :', [node.max_val for node in nodes])
        print('\nindex :', i)
        innodes = [x for x in weight if x[1] == i]
        print('innodes :', innodes)
        print('src-',innodes[0][0], ' : ', nodes[innodes[0][0]].max_val, sep='')
        print('dst-', innodes[0][1], ' : ', nodes[innodes[0][1]].max_val, sep='')

        if len(innodes) < 1:
            continue

        # score = max(prev_node scores + weights)
        print('innode_idx :', [(innode[0], innode[1]) for innode in innodes])
        innodes_val = [nodes[innode[0]].max_val + innode[2] for innode in innodes]

        print('vals :', innodes_val)
        print('max_val of nodes[i] (prev) :', nodes[i].max_val)
        nodes[i].max_val = max(innodes_val)
        nodes[i].innode = innodes[innodes_val.index(nodes[i].max_val)][0]
        print('node_vals :', [node.max_val for node in nodes])
        print('max_val of nodes[i]', nodes[i].max_val)
    
    print('\nnode path :', [str(node.innode) for node in nodes[src+2:]], '\n')
    res = [
        nodes[dst].max_val,
        '->'.join([str(node.innode) for node in nodes[src:]])
    ]

    return res


if __name__ == '__main__':
    # test out InputParser
    params = InputParser(sys.argv[1])

    print(params['source'])
    print(params['sink'])
    print(params['weight'], '\n')

    path_selected = PathSelector(params['source'], params['sink'], params['weight'])
    print(path_selected[0], path_selected[1], sep='\n')
