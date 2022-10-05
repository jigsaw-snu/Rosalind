import sys


'''
    <InputParser>
'''
def InputParser(input_path: str) -> list:
    ret = [0, dict()]

    with open(input_path, 'r') as file:
        # first line refers to number of leaf nodes
        nr_leaf = int(file.readline().rstrip())
        ret[0] = nr_leaf

        while True:
            line = file.readline().rstrip()

            if not line:
                break

            prev_node = int(line[:line.index('-')])
            next_node = int(line[line.index('>')+1:line.index(':')])
            weight = int(line[line.index(':')+1:])

            if prev_node in ret[1].keys():
                ret[1][prev_node].append({next_node:weight})
            else:
                ret[1][prev_node] = [{next_node : weight}]

    return ret


def GetDistance(nr_leaf: int, node_info: dict) -> None:
    for node in node_info.keys():
        near_node = [key for entry in node_info[node] \
                     for (key, value) in entry.items() \
                     if key in range(nr_leaf)]
        print('node:', node, 'near_node:', near_node)


if __name__ == '__main__':
    params = InputParser(sys.argv[1])
    print(params[0], params[1], sep='\n')

    GetDistance(params[0], params[1])
