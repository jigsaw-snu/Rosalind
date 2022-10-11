import time
import sys, copy
from tqdm import tqdm
from collections import deque


'''
    <InputParser>
'''
def InputParser(input_path: str) -> list:
    #ret = [0, dict()]

    # [nr_leaf, {src :[dst]}, {(src, dst) : weight}, [leaf_nodes]]
    ret = [0, dict(), dict(), []]

    with open(input_path, 'r') as file:
        # first line refers to number of leaf nodes
        nr_leaf = int(file.readline().rstrip())
        ret[0] = nr_leaf

        blacklist = []

        while True:
            line = file.readline().rstrip()

            if not line:
                break

            prev_node = int(line[:line.index('-')])
            next_node = int(line[line.index('>')+1:line.index(':')])
            weight = int(line[line.index(':')+1:])

            if prev_node in ret[1].keys():
                ret[1][prev_node].append(next_node)
            else:
                ret[1][prev_node] = [next_node]

            ret[2][(prev_node, next_node)] = weight

            if (prev_node not in ret[3]) and (prev_node not in blacklist):
                ret[3].append(prev_node)
            else:
                # multiple src nodes means not leaf node
                blacklist.append(prev_node)

        # only leaf nodes should be left in ret[3]
        ret[3] = [x for x in ret[3] if x not in blacklist]

    return ret


'''
    <FindPathBFS>
'''
def FindPathBFS(src: int, dst: int, node_info: list) -> list:
    path_history = deque([[src]])  # save all path history in here

    while path_history:
        #print(path_history)
        
        # path_history is a queue, since we want to perform BFS!
        cur_path = path_history.popleft()  # get first path from path_history
        #print(cur_path)
        last_node = cur_path[-1]  # get last node of path
        
        if last_node in cur_path[:-1]:
            continue

        if last_node == dst:  # dst found!
            return cur_path

        try:
            near_nodes = node_info[1][last_node]
        except:
            near_nodes = []

        for near_node in near_nodes:  # spawn new path and add to path_history
            tmp_path = copy.deepcopy(cur_path)
            tmp_path.append(near_node)
            path_history.append(tmp_path)


'''
    <GetDistanceMatrix>
'''
def GetDistanceMatrix(node_info: list) -> None:
    ret_matrix = [[0 for _ in range(node_info[0])] for _ in range(node_info[0])]

    for i in tqdm(range(node_info[0])):
        for j in range(node_info[0]):
            if i == j:
                ret_matrix[i][j] = str(ret_matrix[i][j])
                continue
            
            weight_sum = 0
            node_path = FindPathBFS(node_info[3][i], node_info[3][j], node_info)
            #print('Path', i, 'to', j, ':', node_path)

            for k in range(len(node_path)-1):
                weight_sum += node_info[2][(node_path[k], node_path[k+1])]

            ret_matrix[i][j] = str(weight_sum)
    
    ret_matrix = [' '.join(ret_matrix[i]) for i in range(node_info[0])]

    print('\n', '\n'.join(ret_matrix), sep='')


if __name__ == '__main__':
    params = InputParser(sys.argv[1])
    #print(params[0], params[1], params[2], params[3], '', sep='\n')

    GetDistanceMatrix(params)