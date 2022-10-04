import sys, re
import numpy as np
from copy import deepcopy


'''
    <Function : SubMatParser>
    Parse substitution matrix and save it into objects
'''
def SubMatParser(sub_mat_path: str) -> list:
    ret = [[], []]

    with open(sub_mat_path, 'r') as file:
        while True:
            line = file.readline()  # we'll handle whitespace later

            if not line:
                break

            # A comment starts with '#'
            if '#' in line:
                continue

            # first line is column index
            if not re.findall(r'\d+', line):
                ret[0] += line.split()
                continue

            ret[1].append(list(map(int, line.split()[1:])))

    return [ret[0], np.array(ret[1])]


'''
    <InputParser>
    Parse input file and save into python object
    We consider 2 sequences are given in input file
'''
def InputParser(input_path: str) -> list:
    ret = []

    with open(input_path, 'r') as file:
        ret += [file.readline().rstrip() for _ in range(2)]

    return ret


'''
    <Class : Aligner>
    @variables:
        - seq1, seq2 : target sequences
        - sub_matrix : substitution matrix given as parameter
        - penalty : gap (indel) penalty
        - score_matrix
        - backtrace

    @functions:
        - BacktrackGlobal
        - BacktrackLocal
        - GlobalAligner
        - LocalAligner
'''
class Aligner:
    '''
        <Constructor>
        initialize class member variables

        ** remember using deepcopy to assign same array to different variable
    '''
    def __init__(self, seq1: str, seq2: str, sub_mat_info: list, penalty: int):
        # prepare alignment
        self.seq1, self.seq2 = seq1, seq2
        self.sub_mat_idx = dict((v, k) for (k, v) in enumerate(sub_mat_info[0]))
        self.sub_matrix = sub_mat_info[1]  # np.array
        self.penalty = penalty
        self.score_matrix = \
            np.zeros(shape=(len(self.seq1)+1, len(self.seq2)+1))
        self.backtrace = deepcopy(self.score_matrix)


    '''
        <BacktrackGlobal>
        Global Alignment needs to map sequences from 'end to end'
        You should handle edge cases, unless you'll get incomplete alignment result
    '''
    def BacktrackGlobal(self) -> list:
        ret = [int(self.max_score) ,['', '']]

        i = len(self.seq1)
        j = len(self.seq2)

        while i > 0 and j > 0:
            if self.backtrace[i, j] == 3:  # diagonal movement
                ret[1][0] += self.seq1[i-1]
                ret[1][1] += self.seq2[j-1]
                i -= 1
                j -= 1

            elif self.backtrace[i, j] == 2:  # down movement
                ret[1][0] += self.seq1[i-1]
                ret[1][1] += '-'
                i -= 1

            elif self.backtrace[i, j] == 1:  # right movement
                ret[1][0] += '-'
                ret[1][1] += self.seq2[j-1]
                j -= 1

        # handle edge cases : (i == 0 and j != 0) or (i != 0 and j == 0)
        while i != 0 or j != 0:
            if i == 0 and j != 0:
                ret[1][0] += '-'
                ret[1][1] += self.seq2[j-1]
                j -= 1

            elif i != 0 and j == 0:
                ret[1][0] += self.seq1[i-1]
                ret[1][1] += '-'
                i -= 1

        ret[1][0] = ret[1][0][::-1]
        ret[1][1] = ret[1][1][::-1]

        return ret


    '''
        <BacktrackLocal>
        Local Alignment should start backtracking from max score indices
        You don't need to handle edge cases, since it's not aligning end to end
    '''
    def BacktrackLocal(self, cur_max_score: int, max_score_idx: list) -> list:
        ret = [int(cur_max_score), ['', '']]

        # unlike global alignment, local alignment backtracks from index of max score
        i, j = max_score_idx

        while i > 0 and j > 0 and self.backtrace[i, j] != 0:
            if self.backtrace[i, j] == 3:  # diagonal movement
                ret[1][0] += self.seq1[i-1]
                ret[1][1] += self.seq2[j-1]
                i -= 1
                j -= 1

            elif self.backtrace[i, j] == 2:  # down movement
                ret[1][0] += self.seq1[i-1]
                ret[1][1] += '-'
                i -= 1

            elif self.backtrace[i, j] == 1:  # right movement
                ret[1][0] += '-'
                ret[1][1] += self.seq2[j-1]
                j -= 1

        ret[1][0] = ret[1][0][::-1]
        ret[1][1] = ret[1][1][::-1]

        return ret


    '''
        <GlobalAligner>
        You must initialize first row and column with gap penalty
    '''
    def GlobalAligner(self) -> list:
        # global alignment must initialize first row/col with gap penalty
        self.score_matrix[0, :] = \
            [0 - i*self.penalty for i in range(len(self.seq2)+1)]
        self.score_matrix[:, 0] = \
            [0 - i*self.penalty for i in range(len(self.seq1)+1)]

        # fill out score matrix
        for i in range(1, len(self.seq1)+1):
            for j in range(1, len(self.seq2)+1):
                self.diag = self.score_matrix[i-1, j-1] + \
                                self.sub_matrix[
                                    self.sub_mat_idx[self.seq1[i-1]],
                                    self.sub_mat_idx[self.seq2[j-1]]
                                ]
                self.down = self.score_matrix[i-1, j] - self.penalty
                self.right = self.score_matrix[i, j-1] - self.penalty

                self.max_score = \
                    max([
                        # match or mismatch --> diagonal movement
                        self.diag,
                        # gap in seq2 --> down movement
                        self.down,
                        # gap in seq1 --> right movement
                        self.right
                    ])

                self.score_matrix[i, j] = self.max_score

                if self.max_score == self.diag:
                    self.backtrace[i, j] = 3  # let 3 denotes diagonal movement

                elif self.max_score == self.down:
                    self.backtrace[i, j] = 2  # let 2 denotes down movement

                elif self.max_score == self.right:
                    self.backtrace[i, j] = 1  # let 1 denotes right movement

                # print(self.score_matrix, self.backtrace, '\n', sep='\n')

        return self.BacktrackGlobal()


    '''
        <LocalAligner>
        You don't have to initialize first row and column with gap penalty,
        since it'll automatically select 0 for max value in next cell

        Remember to save max score and max score index
    '''
    def LocalAligner(self) -> list:
        # fill out score matrix
        # also, max score index must be saved
        cur_max_score = 0
        max_score_idx = [len(self.seq1)+1, len(self.seq2)+1]
        for i in range(1, len(self.seq1)+1):
            for j in range(1, len(self.seq2)+1):
                self.diag = self.score_matrix[i-1, j-1] + \
                                self.sub_matrix[
                                    self.sub_mat_idx[self.seq1[i-1]],
                                    self.sub_mat_idx[self.seq2[j-1]]
                                ]
                self.down = self.score_matrix[i-1, j] - self.penalty
                self.right = self.score_matrix[i, j-1] - self.penalty

                self.max_score = \
                    max([
                        # local alignment
                        0,
                        # match or mismatch --> diagonal movement
                        self.diag,
                        # gap in seq2 --> down movement
                        self.down,
                        # gap in seq1 --> down movement
                        self.right
                    ])

                self.score_matrix[i, j] = self.max_score

                if self.max_score > cur_max_score:
                    cur_max_score = self.max_score
                    max_score_idx = [i, j]

                # save backtracking info
                if self.max_score == 0:
                    continue

                elif self.max_score == self.diag:
                    self.backtrace[i, j] = 3  # let 3 denotes diagonal movement

                elif self.max_score == self.down:
                    self.backtrace[i, j] = 2  # let 2 denotes down movement

                elif self.max_score == self.right:
                    self.backtrace[i, j] = 1  # let 1 denotes right movement

                # print(self.score_matrix, self.backtrace, '\n', sep='\n')

        return self.BacktrackLocal(cur_max_score, max_score_idx)


if __name__ == '__main__':
    option = sys.argv[1]  # -l or -g
    sub_mat_info = SubMatParser(sys.argv[2])  # BLOSUM62 or PAM250
    sequences = InputParser(sys.argv[3])

    aligner = Aligner(sequences[0], sequences[1], sub_mat_info, 5)

    if 'g' in option:
        aligned = aligner.GlobalAligner()

    elif 'l' in option:
        aligned = aligner.LocalAligner()

    print(aligned[0], '\n'.join(aligned[1]), sep='\n')
