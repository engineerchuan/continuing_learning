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

exercise_01_02_b()


