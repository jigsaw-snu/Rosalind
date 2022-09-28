import sys, re
import pandas as pd


'''
    <BlosumParser>
    Read Blosum matrix in txt form and save into dataframe form
'''
def BlosumParser(blosum_path: str) -> pd.DataFrame:
    blosum = []
    idx = []
    
    with open(blosum_path, 'r') as blosum:
        while True:
            line = blosum.readline().strip()

            if not line:
                break
            
            if not re.findall(r'\d+', line):
                idx += line
            
            blosum.append(list(map(int, line[1:]))
        
    res = pd.DataFrame(blosum, columns=idx, index=idx)

    return res


'''
    <InputParser>
    Parse input data and save it to dictionary
'''
def InputParser(file_path: str) -> dict:
    params = dict()
    params['seq'] = []

    with open(file_path, 'r') as file:
        while True:
            line = file.readline().rstrip()

            if not line:
                break

            params['seq'].append(line)

    return params

'''
'''
def GlobalAligner(seq1: str, seq2: str, blosum: pd.DataFrame) -> list:
    pass


if __name__ == '__main__':
    blossum = BlosumParser(sys.argv[2])
    params = InputParser(sys.argv[1])
