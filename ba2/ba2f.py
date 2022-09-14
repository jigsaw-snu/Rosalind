import sys
import pandas as pd
import numpy as np
from tqdm import tqdm


'''
    <input_parser>
    parse input file and save data in dictionary
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

            dna += ' '
            dna += line

        params['DNA'] = dna.split()

    return params


'''
    <pssm_builder>
    build profile matrix by given motif sequences
'''
def pssm_builder(motifs: list) -> pd.DataFrame:
    seq_df = pd.DataFrame(list(map(list, motifs)))

    pssm = pd.DataFrame(0.0, index=list('ACGT'), columns=range(len(motifs[0])))

    for i in range(len(motifs[0])):
        for base in list('ACGT'):
            try:
                pssm[i][base] = (seq_df[i].value_counts()[base]+1) / (len(motifs[0])+1)
            except:
                pssm[i][base] = (1.0) / (len(motifs[0])+1)

    return pssm


'''
    <single_motif_finder>
    finds motif from single DNA sequence based on given profile matrix
    helper function for motif_finder
'''
def single_motif_finder(DNA: str, k: int, pssm: pd.DataFrame) -> str:
    kmers = [DNA[i:i+k] for i in range(len(DNA)-k+1)]

    best_score = 0
    motif = kmers[0]
    for kmer in kmers:
        score = 1

        for i in range(k):
            score *= pssm[i][kmer[i]]

        if best_score < score:
            best_score = score
            motif = kmer

    return motif


'''
    <motif_finder>
    find motif in each DNA sequence using profile matrix
'''
def motif_finder(DNA: list, k: int, pssm: pd.DataFrame) -> list:
    motifs = [single_motif_finder(dna, k, pssm) for dna in DNA]

    return motifs


'''
    <score_motifs>
    calculate scores of given motifs using shannon entropy
'''
def score_motifs(motifs: list) -> float:
    # build profile matrix
    pssm = pssm_builder(motifs)

    score = 0
    for i in range(len(motifs[0])):
        score += pssm[i].max() * np.log(pssm[i].max())

    return -score


'''
    <stochastic_motif_search>
    find motif based on randomized algorithm
'''
def stochastic_motif_search(DNA: list, k: int, t: int) -> str:
    # all kmer list of each DNA sequence
    kmers = [[dna[i:i+k] for i in range(len(dna)-k+1)] for dna in DNA]

    # create random index
    random_index = np.random.randint(0, len(DNA[0])-k+1, t).tolist()

    # initialize best_motifs with random index 
    best_motifs = motifs = [kmer[random_index.pop()] for kmer in kmers]

    cnt = 0
    while True:
        cnt += 1

        pssm = pssm_builder(motifs)
        motifs = motif_finder(motifs, k, pssm)

        if cnt == 1000:
            return ('\n').join(best_motifs)
        else:
            best_motifs = motifs
            continue


if __name__ == '__main__':
    params = input_parser(sys.argv[1])
    print(stochastic_motif_search(params['DNA'], params['k'], params['t']))
