import sys, itertools


'''
    <input_parse>
    parse input file and save data into dictionary
'''
def input_parser(file_path: str) -> dict:
    params = dict()

    with open(file_path, 'r') as file:
        # first line refers to k
        line = file.readline().rstrip()
        params['k'] = int(line)

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
    <hamming_distance>
    calculate hamming distance between two sequences
'''
def hamming_distance(seq1: str, seq2: str) -> int:
    if len(seq1) != len(seq2):
        raise Exception("Length of input sequences are different!!")

    score = 0

    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            score += 1

    return score


'''
    <median_string_finder>
    find kmer motif that has minimum hamming distance over every kmers in DNA sequences
'''
def median_string_finder(DNA: list, k: int) -> str:
    kmers = [[dna[i:i+k] for i in range(0, len(dna) - k + 1)] for dna in DNA]

    min_dist = k * len(DNA[0])
    median_motif = ''

    # iterate every possible motifs
    for motif in [''.join(x) for x in itertools.product('ACGT', repeat=k)]:
        dist_sum = 0
        local_motif = ''

        for kmer_list in kmers:
            local_dist = k
            for kmer in kmer_list:
                if hamming_distance(motif, kmer) < local_dist:
                    local_dist = hamming_distance(motif, kmer)
                    local_motif = motif

            dist_sum += local_dist

        if dist_sum < min_dist:
            min_dist = dist_sum
            median_motif = local_motif

    return median_motif


if __name__ == '__main__':
    params = input_parser(sys.argv[1])

    print(median_string_finder(params['DNA'], params['k']))
