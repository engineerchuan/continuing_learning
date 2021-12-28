import sys
sys.path.append("../wk3")
from lib_wk3 import *


def randomized_motif_search(k, t, dnas):
    import random
    n_len = len(dnas[0])
    assert len(dnas) == t
    # initialize the starting motif by random selection
    motifs = []
    for dna in dnas:
        rando_start = random.randint(0, n_len - k)
        motifs.append(dna[rando_start:rando_start+k])
    best_motifs = motifs

    while True:
        profile = compute_profile(motifs, padding=1)
        motifs = []
        for dna in dnas:
            (p, kmer) = most_probable_kmer(dna, k, profile)
            motifs.append(kmer)
        if compute_scores_total(motifs) < compute_scores_total(best_motifs):
            best_motifs = motifs
        else:
            return (compute_scores_total(best_motifs), best_motifs)

def randomized_motif_search_wrapped(k, t, dnas, n_iterations=1):
    best_score_so_far = 1000000000000
    best_motif_so_far = None
    for i in range(n_iterations):
        print(i)
        (best_score, best_motif) = randomized_motif_search(k, t, dnas)
        if best_score < best_score_so_far:
            best_score_so_far = best_score
            best_motif_so_far = best_motif
    return best_motif_so_far

def exercise_01_01_a():
    k = 8
    t = 5
    dnas = ['CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA', 'GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG', 'TAGTACCGAGACCGAAAGAAGTATACAGGCGT', 'TAGATCAAGTTTCAGGTGCACGTCGGTGAACC', 'AATCCACCAGCTCCACGTGCAATGTTGGCCTA']
    with codecs.open('data/dataset_161_5.txt', encoding='utf-8') as fid:
        __ = fid.readline().strip().split(' ')
        k = int(__[0])
        t = int(__[1])
        dnas = fid.readline().strip().split(' ')            
    print(' '.join(randomized_motif_search_wrapped(k, t, dnas, 1000)))

exercise_01_01_a()