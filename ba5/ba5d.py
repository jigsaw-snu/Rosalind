import sys
import numpy as np


'''
    <InputParser>
    Parse input data and save it to a dictionary
'''
def InputParser(file_path: str) -> dict:
    params = dict()

    with open(file_path, 'r') as file:
        # first line indicates source node
        line = file.readline().rstrip()
        params['source'] = int(line)

        # second line indicates sink node
        line = file.readline().rstrip()
        params['sink'] = int(line)

        # from 3rd line, indicates weight information
        # (source)->(sink):(weight)
        # params['weights'] = weighted adjacency matrix
        #
        # **Check**
        # 1. how many nodes are there? --> shape(mat) = (num_nodes) x (num_nodes)
        #    -> we consider number of nodes as size of biggest node
        #    -> check biggest node
        # 1-1. set biggest node = max(source_node, sink_node)
        #      -> source node > sink node might be possible
        # 1-2. while iterating each line, renew biggest node
        #      -> we cannot guarantee source_node or sink_node is the biggest node
        # 2. construct zero matrix with len(biggest_node)
        # 3. substitue each weight to a given position of matrix
        biggest_node = params['sink']
        weight_list = []

        while True:
            line = file.readline().rstrip()

            if not line:
                break

            src, dst, wt = map(int, [
                line[:line.index('-')],
                line[line.index('>')+1:line.index(':')],
                line[line.index(':')+1:]
                    ])

            if biggest_node < src:
                biggest_node = src

            if biggest_node < dst:
                biggest_node = dst

            weight_list.append((src, dst, wt))

        # remember node starts from 0 and biggest_node must be included in the matrix
        params['weights'] = np.zeros(shape=(biggest_node+1, biggest_node+1))
        
        for weight in weight_list:
            src, dst, wt = weight
            
            # we set row as src, col as dst
            params['weights'][src, dst] = wt

    return params
 



if __name__ == '__main__':
    # test out InputParser
    params = InputParser(sys.argv[1])
    print(params['source'])
    print(params['sink'])
    print(params['weights'])
