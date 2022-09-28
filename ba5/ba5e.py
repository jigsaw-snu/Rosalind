import sys
import pandas as pd


''''''
def BlosumParser(blossum_path: str) -> pd.DataFrame:
    blossum = pd.DataFrame()


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



if __name__ == '__main__':
    params = InputParser(sys.argv[1])
