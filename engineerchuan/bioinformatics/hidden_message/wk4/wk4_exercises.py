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


def gibbs_sampler_motif_search(k, t, dnas, N):
    import random
    n_len = len(dnas[0])
    assert len(dnas) == t
    # initialize the starting motif by random selection
    motifs = []
    for dna in dnas:
        rando_start = random.randint(0, n_len - k)
        motifs.append(dna[rando_start:rando_start+k])
    best_motifs = motifs

    for i in range(N):
        # make a random choice and remove a motif
        i_dna_to_remove = random.randint(0, t-1)
        del motifs[i_dna_to_remove]

        # now build the profile
        profile = compute_profile(motifs, padding=1)

        # now roll a biased dice to find the most probable k mer
        dna = dnas[i_dna_to_remove]
        p_s = []
        for i in range(len(dna) - k + 1):
            sliced = dna[i:i+k]
            p = 1.0
            for i in range(len(sliced)):
                p *= profile[sliced[i]][i]
            p_s.append(p)
        
        # roll the biased dice
        i_kmer = random.choices(range(len(dna)-k+1), weights=p_s, k=1)
        #print(i_kmer)
        i_kmer= i_kmer[0]
        #print(motifs)
        #print(dna[i_kmer:(i_kmer+k)])
        motifs.insert(i_dna_to_remove, dna[i_kmer:(i_kmer+k)])
        
        if compute_scores_total(motifs) < compute_scores_total(best_motifs):
            best_motifs = motifs
        #else:
    return (compute_scores_total(best_motifs), best_motifs)

def gibbs_sampler_motif_search_wrapped(k, t, dnas, N, n_wrappers):
    best_score_so_far = 100000000
    best_motif_so_far = None
    for i in range(n_wrappers):
        (score, motif) = gibbs_sampler_motif_search(k, t, dnas, N)
        if score < best_score_so_far:
            best_score_so_far = score
            best_motif_so_far = motif
    return (best_score_so_far, best_motif_so_far)

def exercise_01_02_b():
    k = 8
    t = 5
    N=100
    dnas = ['CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA', 'GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG', 'TAGTACCGAGACCGAAAGAAGTATACAGGCGT', 'TAGATCAAGTTTCAGGTGCACGTCGGTGAACC', 'AATCCACCAGCTCCACGTGCAATGTTGGCCTA']
    with codecs.open('data/dataset_163_4.txt', encoding='utf-8') as fid:
        __ = fid.readline().strip().split(' ')
        k = int(__[0])
        t = int(__[1])
        N = int(__[2])
        dnas = fid.readline().strip().split(' ')            
        print(gibbs_sampler_motif_search_wrapped(k, t, dnas, N, 20))

exercise_01_02_b()
