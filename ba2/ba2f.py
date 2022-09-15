import sys
import pandas as pd
import numpy as np


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
                pssm[i][base] = (seq_df[i].value_counts()[base]+1) / (len(motifs)+4)
            except:
                pssm[i][base] = (1.0) / (len(motifs)+4)

    return pssm


'''
    <single_motif_finder>
    finds motif from single DNA sequence based on given profile matrix
    helper function for motif_finder
    the score used in this function is different from the score used in
    evaluting motifs
'''
def single_motif_finder(dna: str, k: int, pssm: pd.DataFrame) -> str:
    kmers = [dna[i:i+k] for i in range(len(dna)-k+1)]

    best_score = 0
    motif = kmers[0]

    for kmer in kmers:
        score = 1

        for i in range(k):
            score *= pssm[i][kmer[i]]  # score = product of maximum probs

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
    shannon entropy = -sum(sum of P*logP per base for each column)
    lower is better (which refers converging)
'''
def score_motifs(motifs: list) -> float:
    # build profile matrix
    pssm = pssm_builder(motifs)

    score = -sum([sum([pssm[i][base]*np.log2(pssm[i][base]) for base in list('ACGT')]) for i in range(len(motifs[0]))])

    return score


'''
    <stochastic_motif_search>
    find motif based on randomized algorithm
'''
def stochastic_motif_search(DNA: list, k: int, t: int) -> str:
    # all kmer list of each DNA sequence
    total_kmers = [[dna[i:i+k] for i in range(len(dna)-k+1)] for dna in DNA]

    # create random index ; notice that each DNA sequence might has different size
    random_index = [int(np.random.randint(0, len(kmers), 1)) for kmers in total_kmers]
    #random_index = np.random.randint(0, len(DNA[0])-k+1, t).tolist()

    # initialize best_motifs with random index
    best_motifs = curr_motifs = [kmers[random_index.pop()] for kmers in total_kmers[::-1]]

    cnt = 0
    while True:
        cnt += 1
        #print(cnt)

        print('best motifs :', best_motifs)
        pssm = pssm_builder(curr_motifs)
        print('======\nPSSM\n======\n', pssm, sep='')
        curr_motifs = motif_finder(DNA, k, pssm)
        print('='*70)
        print('current motifs :', curr_motifs)
        print(score_motifs(best_motifs), score_motifs(curr_motifs), sep='\t')
        print('\n\n')

        if cnt == 10 or score_motifs(best_motifs) < score_motifs(curr_motifs):
            return ('\n').join(best_motifs)
        else:
            best_motifs = curr_motifs
            continue


if __name__ == '__main__':
    params = input_parser(sys.argv[1])
    print(stochastic_motif_search(params['DNA'], params['k'], params['t']))
