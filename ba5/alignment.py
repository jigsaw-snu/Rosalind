import numpy as np


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
        self.sub_mat_idx = dict((v, k) for (k, v) in enumerate(self.sub_mat_info[0]))
        self.sub_matrix = self.sub_mat_info[1]  # np.array
        self.penalty = penalty
        self.score_matrix = self.backtrace = \
            np.zeros(shape=(len(self.seq1)+1, len(self.seq2)+1))


    def BacktrackGlobal(self, seq1: str, seq2: str, backtrace: np.array) -> list:
        pass


    def BacktrackLocal(self, seq1: str, seq2: str, backtrace: np.array) -> list:
        pass


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

                self.best_score = \
                    max([
                        # match or mismatch --> diagonal movement
                        self.diag,
                        # gap in seq2 --> down movement
                        self.down,
                        # gap in seq1 --> right movement
                        self.right
                    ])

                self.score_matrix[i, j] = self.best_score

                if self.best_score == self.diag:
                    self.backtrace[i, j] = 3  # let 3 denotes diagonal movement

                elif self.best_score == self.down:
                    self.backtrace[i, j] = 2  # let 2 denotes down movement

                elif self.best_score == self.right:
                    self.backtrace[i, j] = 1  # let 1 denotes right movement

        return self.BacktrackGlobal(self.seq1, self.seq2, self.backtrace)


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

                self.best_score = \
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

                if self.best_score == 0:
                    continue

                elif self.best_score == self.diag:
                    self.backtrace[i, j] = 3  # let 3 denotes diagonal movement

                elif self.best_score == self.down:
                    self.backtrace[i, j] = 2  # let 2 denotes down movement

                elif self.best_score == self.right:
                    self.backtrace[i, j] = 1  # let 1 denotes right movement

        return self.BacktrackLocal(self.seq1, self.seq2, self.backtrace)


if __name__ == '__main__':
    pass
