from os import error
import sys
import pandas as pd
import numpy as np
from tqdm import tqdm


'''
    <input_parser>
    parse input file and save into dictionary
'''
def input_parser(file_path: str) -> dict:
    params = dict()

    with open(file_path, 'r') as file:
        # first line refers to k and t
        line = file.readline().rstrip()
        params['k'], params['t'] = list(map(int, line.split()))

        # from second line, it refers to multiple DNA sequences
        dna = ''

        while True:
            line = file.readline().rstrip()

            if not line:
                break

            # for safety, we rather not to add lines to list directly
            dna += ' '
            dna += line

        params['DNA'] = dna.split()

    return params


'''
    <pssm_builder>
    build profile matrix with given DNA sequences
    this time, add pseudocounts to avoid 0 values
'''
def pssm_builder(DNA: list, k: int):
    seq_df = pd.DataFrame(list(map(list, DNA)))
    pssm = pd.DataFrame(0.0, index=list('ACGT'), columns=range(k))

    for i in range(k):
        for base in list('ACGT'):
            try:
                pssm[i][base] = (seq_df[i].value_counts()[base]+1) / (len(DNA)+4)
            except:
                pssm[i][base] = (1.0) / (len(DNA)+4)

    return pssm


'''
    <motif_profiler>
    find motif from given DNA sequence using profile matrix
'''
def motif_profiler(DNA: str, k: int, pssm: pd.DataFrame) -> str:
    kmers = [DNA[i:i+k] for i in range(len(DNA)-k+1)]

    # if all the scores are same, then just return first kmer
    best_motif = kmers[0]
    best_score = 0

    for kmer in kmers:
        score = 1

        for i in range(k):
            score *= pssm[i][kmer[i]]

        if score > best_score:
            best_score = score
            best_motif = kmer

    return best_motif


'''
    <score_motifs>
    calculate shannon entropy score of given set of motifs
    ** greater is better **
'''
def score_motifs(motifs: list):
    # build probability matrix, which is same as N-mer pssm
    probs = pssm_builder(motifs, len(motifs[0]))

    entropy = 0
    for i in range(len(motifs[0])):
        entropy += probs[i].max() * np.log2(probs[i].max())

    return -entropy


'''
    <greedy_motif_finder>
    find motif using greedy algorithm based on profile matrix
    this time, we add pseudocounts to profile matrix
'''
def greedy_motif_finder(DNA: list, k: int, t: int) -> str:
    kmers = [[dna[i:i+k] for i in range(len(dna)-k+1)] for dna in DNA]

    # initialize best_motifs with first kmer of each DNA sequence
    best_motifs = [kmer[0] for kmer in kmers]

    # iterate over all kmers in first DNA sequence
    for kmer in tqdm(kmers[0]):
        motifs = [kmer]

        for i in range(1, t):
            pssm = pssm_builder(motifs[:i], k)
            new_motif = motif_profiler(DNA[i], k, pssm)
            motifs.append(new_motif)

        # compare score between best_motifs and current_motifs
        if score_motifs(best_motifs) < score_motifs(motifs):
            best_motifs = motifs

    return ('\n').join(best_motifs)


if __name__ == '__main__':
    params = input_parser(sys.argv[1])

    print(greedy_motif_finder(params['DNA'], params['k'], params['t']))
