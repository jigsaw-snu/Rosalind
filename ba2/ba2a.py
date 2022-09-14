import sys


'''
    <input_parser>
    parse input file and save data into dictionary
'''
def input_parser(file_path: str) -> dict:
    params = dict()

    with open(file_path, 'r') as file:
        # first line refers to k and d
        line = file.readline().rstrip()
        params['k'], params['d'] = list(map(int, line.split()))

        # from second line, it refer to DNA sequences
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
    <HammingDistance>
    Calculate Hamming Distance of Two sequences
    Length of input sequences must be equal
'''
def hamming_distance(seq1: str, seq2: str) -> int:
    if len(seq1) != len(seq2):
        raise Exception("Length of input sequences are different!")

    score = 0

    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            score += 1

    return score


'''
    <neighbors_Dd>
    Find sequence that has hamming distance same or less than d.
'''
def neighbors_Dd(seq: str, d: int) -> list:
    neighbors = []

    if d == 0:
        return [seq]

    if len(seq) < 2:  # base case
        return ['A', 'C', 'G', 'T']

    init_base = seq[0]
    suffix = seq[1:]

    suffix_neighbors = neighbors_Dd(suffix, d)  # recursive case

    for suff_nb in suffix_neighbors:
        if hamming_distance(suffix, suff_nb) < d:
            neighbors += [x + suff_nb for x in ['A', 'C', 'G', 'T']]
        else:
            neighbors += [init_base + suff_nb]

    return neighbors


'''
    <motif_enumeration>
    A Brute Force algorithm to find motifs
'''
def motif_enumeration(DNA: list, k: int, d: int) -> str:
    motifs = []

    kmers = [[dna[i:i+k] for i in range(0, len(dna) - k +1)] for dna in DNA]

    for kmer in kmers[0]:
        k_neighbors = neighbors_Dd(kmer, d)

        for kn in k_neighbors:
            dna_token = 0

            for other_kmers in kmers[1:]:
                kmer_token = 0

                for other_kmer in other_kmers:
                    if hamming_distance(kn, other_kmer) > d:
                        continue
                    else:  # hamming distnace < d or == d is acceptable
                        kmer_token += 1

                if kmer_token > 0:
                    dna_token += 1

            if dna_token == len(DNA) - 1:  # motif should appear in each sequences
                motifs.append(kn)

    motifs = list(set(motifs))

    return ' '.join(motifs)


if __name__ == '__main__':
    params = input_parser(sys.argv[1])
    print(motif_enumeration(params['DNA'], params['k'], params['d']))
