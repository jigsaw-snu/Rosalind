import sys, re
import numpy as np


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
    <Class : Aligner>
    @variables:
        - seq1, seq2 : target sequences
        - sub_matrix : substitution matrix given as parameter
        - penalty : gap (indel) penalty
        - score_matrix
        - backtrace

    @functions:
        - GlobalAligner
        - LocalAligner
'''
class Aligner:
    def __init__(self, seq1: str, seq2: str, sub_mat_info: list, penalty: int):
        # prepare alignment
        self.seq1, self.seq2 = seq1, seq2
        self.sub_mat_idx = dict((v, k) for (k, v) in enumerate(sub_mat_info[0]))
        self.sub_matrix = sub_mat_info[1]  # np.array
        self.penalty = penalty
        self.score_matrix = self.backtrace = \
            np.zeros(shape=(len(self.seq1)+1, len(self.seq2)+1))


    def BacktrackGlobal(self) -> list:
        ret = [self.max_score ,['', '']]

        i = len(self.seq1)
        j = len(self.seq2)

        while i > 0 and j > 0:
            if self.backtrace[i, j] == 3:  # diagonal movement
                ret[1][0] += self.seq1[i-1]
                ret[1][1] += self.seq2[j-1]
                i -= 1
                j -= 1

            if self.backtrace[i, j] == 2:  # down movement
                ret[1][0] += self.seq1[i-1]
                ret[1][1] += '-'
                i -= 1

            if self.backtrace[i, j] == 1:  # right movement
                ret[1][0] += '-'
                ret[1][1] += self.seq2[j-1]
                j -= 1

        ret[1][0] = ret[1][0][::-1]
        ret[1][1] = ret[1][1][::-1]

        return ret


    def BacktrackLocal(self) -> list:
        ret = [self.max_score, []]

        i = len(self.seq1)
        j = len(self.seq2)

        while i > 0 and j > 0 and self.backtrace[i, j] != 0:
            if self.backtrace[i, j] == 3:  # diagonal movement
                ret[1][0] += self.seq1[i-1]
                ret[1][1] += self.seq2[j-1]
                i -= 1
                j -= 1

            if self.backtrace[i, j] == 2:  # down movement
                ret[1][0] += self.seq1[i-1]
                ret[1][1] += '-'
                i -= 1

            if self.backtrace[i, j] == 1:  # right movement
                ret[1][0] += '-'
                ret[1][1] += self.seq2[j-1]
                j -= 1

        ret[1][0] = ret[1][0][::-1]
        ret[1][1] = ret[1][1][::-1]

        return ret


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

                print(self.diag, self.down, self.right)

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

        print(self.score_matrix, '\n')

        return self.BacktrackGlobal()


    def LocalAligner(self) -> list:
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
                        # local alignment
                        0,
                        # match or mismatch --> diagonal movement
                        self.diag,
                        # gap in seq2 --> down movement
                        self.down,
                        # gap in seq1 --> down movement
                        self.right
                    ])

                if self.max_score == 0:
                    continue

                elif self.max_score == self.diag:
                    self.backtrace[i, j] = 3  # let 3 denotes diagonal movement

                elif self.max_score == self.down:
                    self.backtrace[i, j] = 2  # let 2 denotes down movement

                elif self.max_score == self.right:
                    self.backtrace[i, j] = 1  # let 1 denotes right movement

        print(self.score_matrix, '\n')

        return self.BacktrackLocal()


if __name__ == '__main__':
    sub_mat_info_gl = SubMatParser(sys.argv[1])  # BLOSUM62
    sub_mat_info_lc = SubMatParser(sys.argv[2])  # PAM250
    print(sub_mat_info_gl[0], sub_mat_info_gl[1], sep='\n')
    print('\n')
    print(sub_mat_info_lc[0], sub_mat_info_lc[1], sep='\n')
    print('\n')

    # Test Global Alignment
    aligner_gl = Aligner('PLEASANTLY', 'MEANLY', sub_mat_info_gl, 5)
    aligned_gl = aligner_gl.GlobalAligner()
    print(aligned_gl[0], '\n'.join(aligned_gl[1]), sep='\n')
    print('\n')

'''
    # Test Local Alignment
    aligner_lc = Aligner('MEANLY', 'PENALTY', sub_mat_info_lc, 5)
    aligned_lc = aligner_lc.LocalAligner()
    print(aligned_lc[0], '\n'.join(aligned_lc[1]), sep='\n')
'''
