import re, sys


class Node:
    def __init__(self, node_idx):
        self.max_val = 0
        self.innode = [str(node_idx)]


'''
    <StepikParser>
    Parse input data from stepilk.org
'''
def StepikParser(file_path: str) -> dict:
    params = dict()

    with open(file_path, 'r') as file:
        # first line implies source and sink node
        line = file.readline().rstrip()
        params['source'], params['sink'] = list(map(int, line.split()))

        # from second line, weight information exist
        params['weight'] = []

        while True:
            line = file.readline().rstrip()

            if not line:
                break

            src, dst, wt = list(map(int, line.split()))
            params['weight'].append((src, dst, wt))

    return params


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
        params['weight'] = []

        while True:
            line = file.readline().rstrip()

            if not line:
                break

            src, dst, wt = list(map(int, re.findall(r'\d+', line)))
            params['weight'].append((src, dst, wt))

    return params


'''
    <PathSelector>
    Select Best Path from source node to sink node
'''
def PathSelector(src: int, dst: int, weight: list) -> list:
    # this solution only works when weight list is sorted by its source node
    weight.sort()

    # add different Node objects to nodes list
    # make list index same as iterator below --> make (dst+1) nodes
    nodes = [Node(x) for x in range(dst+1)]

    for i in range(src+1, dst+1):
        # (src, dst, wt)
        innodes = [x for x in weight if x[1] == i]

        if len(innodes) < 1:
            continue

        # score = max(prev_node scores + weights)
        innodes_score = [
            nodes[innode[0]].max_val + innode[2] \
            if str(src) in nodes[innode[0]].innode or src == innode[0] else 0 \
            for innode in innodes
        ]

        max_score = max(innodes_score)

        if max_score == 0:
            continue

        nodes[i].max_val = max(innodes_score)
        max_score_pair = innodes[innodes_score.index(max(innodes_score))]

        nodes[i].innode = nodes[max_score_pair[0]].innode + nodes[i].innode

    res = [nodes[dst].max_val, '->'.join(nodes[dst].innode)]

    return res


if __name__ == '__main__':
    # test out InputParser
    if sys.argv[1] == '-s':
        params = StepikParser(sys.argv[2])
    else:
        params = InputParser(sys.argv[1])

    path_selected = PathSelector(params['source'], params['sink'], params['weight'])
    print(path_selected[0], path_selected[1], sep='\n')
