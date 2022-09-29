import sys, re
import numpy as np


'''
    <PAMParser>
    Parse txt format of PAM matrix
'''
def PAMParser(pam_path: str) -> list:
    contents = []
    p_idx = []

    with open(pam_path, 'r') as pam:
        while True:
            line = pam.readline()
            
            if not line:
                break
            
            # line starts with '#' is a comment
            if '#' in line:
                continue
            
            if not re.findall(r'\d+', line):
                p_idx += line.split()
                continue
            
            
            contents.append(list(map(int, line.split()[1:])))

    contents = np.array(contents)
    
    return (p_idx, contents)


'''
    <InputParser>
    Parse input data and save it into a dictionary
'''
def InputParser(file_path: str) -> dict:
    params = dict()
    
    with open(file_path, 'r') as file:
        params['sequence'] = list(map(lambda x: x.rstrip(), file.readlines()))

    return params


'''
    <LocalAligner>
    Alingn two given sequences
    Find the subsequence that maximizes the score
'''
def LocalAligner(sequence: list, pam: list, penalty: int) -> list:
    # Srttings for PAM matrix
    p_idx = dict((v, k) for v, k in enumerate(pam[0]))
    p_mat = pam[1]
    
    # Build Score Matrix & Backtrack Matrix
    # you must include padding row & column (first row/col)
    score_matrix = np.zeros(shape=(len(sequence[0])+1, len(sequence[1])+1))

    # Initialize first row/column with gap or mismatch penalty 
    score_matrix[:, 0] = [0 - i*penalty for i in range(len(sequence[0])+1)]
    score_matrix[0, :] = [0 - i*penalty for i in range(len(sequence[1])+1)]

    for i in range(1, len(sequence[0])+1):
        for j in range(1, len(sequence[1])+1):
            score_matrix[i, j] = max([
                # quick path for local alignment
                0,
            # diagonal movement (match or mismatch)
            score_matrix[i-1, j-1] + p_mat[p_idx[sequence[0][i-1]], p_idx[sequence[1][j-1]]],
            
            ])


if __name__ == '__main__':
    #p_idx, pam = PAMParser(sys.argv[2])
    #print(p_idx, pam, sep='\n')
    
    params = InputParser(sys.argv[1])
    print(params['sequence'])
    