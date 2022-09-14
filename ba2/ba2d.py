import sys
import pandas as pd
from tqdm import tqdm


'''
    <input_parser>
    parse input file and save data into a dictionary
'''
def input_parser(file_path: str) -> dict:
    params = dict()

    with open(file_path, 'r') as file:
        # first line refers to k and t
        line = file.readline().rstrip()
        params['k'], params['t'] = list(map(int, line.split()))

        # from second line, it refers to DNA sequences
        dna = ''
        while True:
            line = file.readline().rstrip()

            if not line:
                break

            dna += ' '
            dna += line

        params['DNA'] = dna.split()  # to handle case when dna sequences are given as single string

    return params


'''
    <pssm_builder>
    build proflie matrix with given sequences
'''
def pssm_builder(DNA: list, k: int, t: int) -> pd.DataFrame:
    pssm = pd.DataFrame(0.0, index=['A', 'C', 'G', 'T'], columns=list(range(k)))

    dna_df = pd.DataFrame(list(map(list, DNA)))

    for i in range(k):
        for base in ['A', 'C', 'G', 'T']:
            try:
                pssm[i][base] = dna_df[i].value_counts()[base] / t
            except:
                pssm[i][base] = 0.0

    return pssm


'''
    <motif_profiler>
    find motif of given sequence based on profile matrix
'''
def motif_profiler(DNA: str, k: int, pssm: pd.DataFrame) -> str:
    kmers = [DNA[i:i+k] for i in range(0, len(DNA) - k + 1)]

    global_score = 0
    motif = kmers[0]

    for kmer in kmers:
        local_score = 1

        for i in range(k):
            local_score *= pssm[i][kmer[i]]

        if local_score > global_score:
            global_score = local_score
            motif = kmer

    return motif


'''
    <score_motifs>
    calculate count based scoring of given motifs
    ** lower is better **
'''
def score_motifs(motifs: pd.DataFrame) -> int:
    score_sum = 0

    for i in range(len(motifs.columns)):
        score_sum += max(motifs[i].value_counts())

    return len(motifs.columns) * len(motifs.columns) - score_sum


'''
    <greedy_motif_finder>
    find motif using greedy algorithm with profile matrix
'''
def greedy_motif_finder(DNA: list, k: int, t: int) -> str:
    kmers = [[dna[i:i+k] for i in range(0, len(DNA[0]) - k + 1)] for dna in DNA]

    best_motifs = [kmer_list[0] for kmer_list in kmers]

    for kmer in tqdm(kmers[0]):
        motifs = [kmer]

        for j in range(1, t):
            pssm = pssm_builder(motifs[:j], k, j)
            motif_j = motif_profiler(DNA[j], k, pssm)

            if motif_j == '':
                break

            motifs.append(motif_j)

        if len(motifs) < t:
            continue

        if score_motifs(pd.DataFrame(list(map(list, motifs)))) < score_motifs(pd.DataFrame(list(map(list, best_motifs)))):
            best_motifs = motifs

    return '\n'.join(best_motifs)


if __name__ == '__main__':
    params = input_parser(sys.argv[1])
    print(greedy_motif_finder(params['DNA'], params['k'], params['t']))
