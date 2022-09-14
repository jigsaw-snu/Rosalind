import sys
import pandas as pd


'''
    <input parser>
    parse input file and save data into a dictionary
'''
def input_parser(file_path: str) -> dict:
    params = dict()

    with open(file_path, 'r') as file:
        # first line refers to a single DNA sequence
        line = file.readline().rstrip()
        params['DNA'] = line

        # second line refers to k
        line = file.readline().rstrip()
        params['k'] = int(line)

        # from third line, it refers to profile matrix
        lines = []
        while True:
            line = file.readline().rstrip()

            if not line:
                break

            lines.append(list(map(float, line.split())))

        params['PSSM'] = pd.DataFrame(lines)
        params['PSSM'].columns = list(range(params['k']))
        params['PSSM'].index = ['A', 'C', 'G', 'T']

    return params


'''
    <motif_profiler>
    find motif based on profile matrix
'''
def motif_profiler(DNA: str, k: int, pssm: pd.DataFrame) -> str:
    kmers = [DNA[i:i+k] for i in range(0, len(DNA) - k + 1)]

    global_score = 0
    motif = ''

    for kmer in kmers:
        local_score = 1
        for i in range(k):
            local_score *= pssm[i][kmer[i]]

        if local_score > global_score:
            global_score = local_score
            motif = kmer

    return motif


if __name__ == '__main__':
    params = input_parser(sys.argv[1])
    print(motif_profiler(params['DNA'], params['k'], params['PSSM']))
