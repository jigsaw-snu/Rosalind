import sys
import numpy as np


'''
    <InputParser>
    Parse input data
'''
def InputParser(file_path: str) -> dict:
    params = dict()

    with open(file_path, 'r') as file:
        # first line implies match reward and mismatch / inel penalty
        line = file.readline().rstrip()

        params['match'], params['mismatch'], params['indel'] = list(map(int, line.split()))

        # second and third line are nucleotide strings
        params['sequences'] = []

        for _ in range(2):
            params['sequences'].append(file.readline().rstrip())

    return params


'''
    <SeqAligner>
    Global Alignment with match reward and mismatch / indel penalty
'''
def SeqAligner(seq1: str, seq2: str, \
               match_reward: int, mismatch_penalty: int, \
               indel_penalty: int) -> list:

    

    for i in range(len(seq1):
        for j in range(len(seq2)):
            pass
