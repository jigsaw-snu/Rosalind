import sys, re, time
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
    p_idx = dict((v, k) for k, v in enumerate(pam[0]))
    p_mat = pam[1]

    # Build Score Matrix & Backtrack Matrix
    # you must include padding row & column (first row/col)
    score_matrix = np.zeros(shape=(len(sequence[0])+1, len(sequence[1])+1))

    # Initialize first row/column with gap or mismatch penalty 
    #score_matrix[:, 0] = [0 - i*penalty for i in range(len(sequence[0])+1)]
    #score_matrix[0, :] = [0 - i*penalty for i in range(len(sequence[1])+1)]

    # Build Backtrace matrix    ** shape is same as score_matrix 
    backtrace = np.zeros(shape=(len(sequence[0])+1, len(sequence[1])+1))

    # Fill out the score_matrix and backtrace matrix
    #
    # remember row corresponds to seqeunce[0] and
    # column corresponds to sequence[1]

    max_init = 0  # initializaing maximum value
    max_loc = (0, 0)  # initializing location of maximum value

    for i in range(1, len(sequence[0])+1):
        for j in range(1, len(sequence[1])+1):
            # diagonal movement (match or mismatch) --> prev_score + PAM_score (match / mismatch)
            diag = score_matrix[i-1, j-1] + \
                p_mat[p_idx[sequence[0][i-1]], p_idx[sequence[1][j-1]]]

            # right movement (gap in sequence[0]) --> prev_score + indel_penalty
            right = score_matrix[i, j-1] - penalty

            # down movement (gap in sequence[1]) --> prev_score + indel_penalty
            down = score_matrix[i-1, j] - penalty

            # 0 is crucial in local alignment
            max_score = max([0, diag, right, down])
            score_matrix[i, j] = max_score

            if max_score > max_init:
                max_init = max_score
                max_loc = (i, j)

            if max_score == 0:
                continue

            elif max_score == diag:
                backtrace[i, j] = 3  # let 3 as diagonal movement

            elif max_score == right:
                backtrace[i, j] = 2  # let 2 as right movement

            elif max_score == down:
                backtrace[i, j] = 1  # let 1 as down movement

    #print(score_matrix, '\n')
    #print(backtrace, '\n')
    #print(max_init, max_loc, '\n')

    # Backtracing
    res = ['', '', int(max_init)]

    i, j = max_loc

    while i >= 0 and j >= 0 and backtrace[i, j] != 0:
        #print('i:', i, 'j:', j, 'backtrace:', backtrace[i, j])
        #time.sleep(0.5)
        if backtrace[i, j] == 3:  # diagonal movement : match or mismatch
            res[0] += sequence[0][i-1]
            res[1] += sequence[1][j-1]
            i -= 1
            j -= 1

        elif backtrace[i, j] == 2:  # right movement : gap in sequence[0]
            res[0] += '-'
            res[1] += sequence[1][j-1]
            j -= 1

        elif backtrace[i, j] == 1:  # down movement : gap in sequence[1]
            res[0] += sequence[0][i-1]
            res[1] += '-'
            i -= 1

    print(backtrace)

    return res


if __name__ == '__main__':
    p_idx, pam = PAMParser(sys.argv[2])
    #print(p_idx, pam, sep='\n')

    params = InputParser(sys.argv[1])
    #print(params['sequence'], '\n')

    aligned = LocalAligner(params['sequence'], [p_idx, pam], 5)
    print(aligned[2], aligned[0][::-1], aligned[1][::-1], sep='\n')
