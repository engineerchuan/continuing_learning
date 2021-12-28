import codecs

def longest_common_substring(string1, string2):
    # initiate the matrix, first index will be for string1
    # second index for string 2
    # each matrix element will be the longest length
    # and the index in string1
    mat = [[(0, "") for i in range(len(string2)+1)] for j in range(len(string1)+1)]
    for i in range(len(string1)):
        for j in range(len(string2)):
            if string1[i] == string2[j]:
                (max_diagonal, string_so_far) = mat[i][j]
                mat[i+1][j+1] = (max_diagonal+1, string_so_far + string1[i])
            else:
                mat[i+1][j+1] = (0, "")
    # find the max
    max_so_far = 0
    lcs_so_far = ""
    for i in range(len(mat)):
        for j in range(len(mat[0])):
            (max_len, lcs) = mat[i][j]
            if max_len > max_so_far:
                max_so_far = max_len
                lcs_so_far = lcs
    return (max_so_far, lcs_so_far)

def exercise_01_02_a():
    a1 = 'atgaccgggatactgataaaaaaaagggggggggcgtacacattagataaacgtatgaagtacgttagactcggcgccgccg'
    a2 = 'acccctattttttgagcagatttagtgacctggaaaaaaaatttgagtacaaaacttttccgaataaaaaaaaaggggggga'
    a3 = 'tgagtatccctgggatgacttaaaaaaaagggggggtgctctcccgatttttgaatatgtaggatcattcgccagggtccga'
    (len_1_2, lcs) = longest_common_substring(a1, a2)
    print(longest_common_substring(lcs, a3))

def neighbors(pattern, d):
    if d == 0:
        return [pattern]
    else:
        recursive = neighbors(pattern, d-1)
        ans = [pattern]
        for patt in recursive:
            for i in range(len(patt)):
                s = [x for x in ['A','T','G','C'] if x != patt[i]]
                for letter in s:
                    newpatt = [x for x in patt]
                    newpatt[i] = letter
                    ans.append(''.join(newpatt))
        return list(set(ans))

def hamming_distance(genome1, genome2):
    if len(genome1) != len(genome2):
        raise ValueError('the two inputs must be same length')
    mismatches = 0
    for i in range(len(genome1)):
        if genome1[i] != genome2[i]:
            mismatches += 1
    return mismatches

def motif_find_brute_force(k, d, dnas):
    # collect all unique starter patterns
    # O(k * len(dnas))?
    patterns = []
    for dna in dnas:
        for i in range(len(dna) - k + 1):
            patterns.append(dna[i:i+k])
    patterns = list(set(patterns))

    # then collect the neighbors

    expanded = [neighbors(pattern, d) for pattern in patterns]
    flattened_expanded = list(set([item for sublist in expanded for item in sublist]))

    final_patterns = []
    for pattern in flattened_expanded:
        all_match_within = True
        for dna in dnas:
            dna_match = False
            for i in range(len(dna) - len(pattern) + 1):
                sliced = dna[i:i+len(pattern)]
                if hamming_distance(sliced, pattern) <= d: # consider this okay
                    dna_match = True
                    break
            if not dna_match:
                all_match_within = False
                break
        if all_match_within:
            final_patterns.append(pattern)
    return final_patterns

def compute_scores(motif):
    # motif is a list of lists
    # dimension 1 is # of motifs
    # dimension 2 is length of motifs
    # return a score of form A,T,C,G as a list of lists
    n_strands = len(motif)
    n_len = len(motif[0])
    scores = {
        'A': [0 for i in range(n_len)],
        'C': [0 for i in range(n_len)],
        'G': [0 for i in range(n_len)],
        'T': [0 for i in range(n_len)],
    }
    for i_motif in range(n_strands):
        for j in range(n_len):
            scores[motif[i_motif][j]][j] += 1
    return scores

def compute_scores_total(motif):
    score = 0
    scores = compute_scores(motif)
    for i in range(len(scores['A'])):
        t = [scores['A'][i], scores['C'][i],scores['T'][i],scores['G'][i]]
        score += sum(t) - max(t)
    return score


def compute_profile(motif, padding=0):
    # motif is a list of lists
    # dimension 1 is # of motifs
    # dimension 2 is length of motifs
    # return a score of form A,T,C,G as a list of lists
    scores = compute_scores(motif)
    n_len = len(scores['A'])
    profile = {'A':[], 'G':[], 'C':[], 'T':[]}
    for i in range(n_len):
        total = scores['A'][i] + scores['T'][i] + scores['C'][i] + scores['G'][i]
        for v in ['A', 'T', 'C', 'G']:
            profile[v].append((padding + scores[v][i]) *1.0 / (total + 4 * padding))
            
    return profile


import math

def compute_entropy(motif):
    scores = compute_scores(motif)
    n_len = len(scores['A'])
    total_entropy = 0
    for i in range(n_len):
        total = scores['A'][i] + scores['T'][i] + scores['C'][i] + scores['G'][i]
        for v in ['A', 'T', 'C', 'G']:
            p = scores[v][i] *1.0 / total
            if p > 0:
                total_entropy += p * math.log(p, 2)
    total_entropy = total_entropy * -1
    return total_entropy

def kmerdist(pattern, dnas):
    def helper(pattern, dna):
        minDist = 10000000
        for i in range(len(dna) - len(pattern) + 1):
            computed = hamming_distance(pattern, dna[i:i+len(pattern)])
            if computed < minDist:
                minDist = computed
        return minDist
    return sum([helper(pattern, dna) for dna in dnas])


def median_string(k, dnas):
    # generate all possible k-mers from dnas

    patterns = []
    for dna in dnas:
        for i in range(len(dna) - k + 1):
            patterns.append(dna[i:i+k])
    patterns = list(set(patterns))
    distance = 10000000000000
    median_pattern = None
    for pattern in patterns:
        computed_dist = kmerdist(pattern, dnas)
        if distance > computed_dist:
            distance = computed_dist
            median_pattern = pattern
    return (distance, median_pattern)







def exercise_01_02_b():
    # motif finding
    # Input: Integers k and d, followed by a space-separated collection of strings Dna.
    # Output: All (k, d)-motifs in Dna.
    with codecs.open('data/dataset_156_8.txt', encoding='utf-8') as fid:
        __ = fid.readline().strip().split(' ')
        k = int(__[0])
        d = int(__[1])
        dnas = fid.readline().strip().split(' ')
        print(' '.join(motif_find_brute_force(k, d, dnas)))



def exercise_01_03_c():
    motif = [
        'TCGGGGGTTTTT',
        'CCGGTGACTTAC',
        'ACGGGGATTTTC',
        'TTGGGGACTTTT',
        'AAGGGGACTTCC',
        'TTGGGGACTTCC',
        'TCGGGGATTCAT',
        'TCGGGGATTCCT',
        'TAGGGGAACTAC',
        'TCGGGTATAACC',        
    ]
    print(compute_scores(motif))
    print(compute_entropy(motif))

def exercise_01_04_a():
    k = 3
    dnas = ['AAATTGACGCAT','GACGACCACGTT','CGTCAGCGCCTG','GCTGAGCACCGG','AGTTCGGGACAG']

    with codecs.open('data/dataset_158_9.txt', encoding='utf-8') as fid:
        k = int(fid.readline().strip())
        dnas = fid.readline().strip().split(' ')
        print(median_string(k, dnas))

def most_probable_kmer(text, k, profile):
    # profile is a map.
    maxp = -1
    found = None
    for i in range(len(text) - k + 1):
        sliced = text[i:i+k]
        p = 1.0
        for i in range(len(sliced)):
            p *= profile[sliced[i]][i]
        if p > maxp:
            maxp = p
            found = sliced

    return (p, found)

def exercise_01_05_a():
    test_profile = {
        'A': [0.2, 0.2, 0.3, 0.2, 0.3],
        'C': [0.4, 0.3, 0.1, 0.5, 0.1],
        'G': [0.3, 0.3, 0.5, 0.2, 0.4],
        'T': [0.1, 0.2, 0.1, 0.1, 0.2],
    }

    with codecs.open('data/dataset_159_3.txt', encoding='utf-8') as fid:
        
        text = fid.readline().strip()
        k = int(fid.readline().strip())
        test_profile = dict()
        test_profile['A'] = [float(x) for x in fid.readline().strip().split(' ')]
        test_profile['C'] = [float(x) for x in fid.readline().strip().split(' ')]
        test_profile['G'] = [float(x) for x in fid.readline().strip().split(' ')]
        test_profile['T'] = [float(x) for x in fid.readline().strip().split(' ')]
        print(most_probable_kmer(text, k, test_profile))

def greedy_motif_search(dnas, k, t, padding=0):
    # initialize the motifs
    best_motifs = [dna[:k] for dna in dnas]
    for i in range(len(dnas[0]) - k + 1):
        motifs = [dnas[0][i:i+k]] # seed it with the first one
        for j in range(1, t):
            #print(motifs)
            profile = compute_profile(motifs, padding=padding)
            (p, kmer) = most_probable_kmer(dnas[j], k, profile)
            motifs.append(kmer)
        #print(motifs)
        if compute_scores_total(motifs) < compute_scores_total(best_motifs):
            best_motifs = motifs
    return best_motifs

def exercise_01_05_b():

    k = 3
    t = 5
    dnas = ['GGCGTTCAGGCA', 'AAGAATCAGTCA', 'CAAGGAGTTCGC', 'CACGTCAATCAC', 'CAATAATATTCG']

    with codecs.open('data/dataset_159_5.txt', encoding='utf-8') as fid:
        __ = fid.readline().strip().split(' ')
        k = int(__[0])
        t = int(__[1])
        dnas = fid.readline().strip().split(' ')        
        print(' '.join(greedy_motif_search(dnas, k, t)))

#exercise_01_05_b()


def exercise_01_05_c():

    k = 3
    t = 5
    dnas = ['GGCGTTCAGGCA', 'AAGAATCAGTCA', 'CAAGGAGTTCGC', 'CACGTCAATCAC', 'CAATAATATTCG']

    with codecs.open('data/dataset_160_9.txt', encoding='utf-8') as fid:
        __ = fid.readline().strip().split(' ')
        k = int(__[0])
        t = int(__[1])
        dnas = fid.readline().strip().split(' ')        
        print(' '.join(greedy_motif_search(dnas, k, t, padding=1)))

def quiz():
    print(median_string(7, ['CTCGATGAGTAGGAAAGTAGTTTCACTGGGCGAACCACCCCGGCGCTAATCCTAGTGCCC', 'GCAATCCTACCCGAGGCCACATATCAGTAGGAACTAGAACCACCACGGGTGGCTAGTTTC', 'GGTGTTGAACCACGGGGTTAGTTTCATCTATTGTAGGAATCGGCTTCAAATCCTACACAG']))
    print(0.4 * 0.3 * 1.0 * 0.4 * 0.5 * 0.1)