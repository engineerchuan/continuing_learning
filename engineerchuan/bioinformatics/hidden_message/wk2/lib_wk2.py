def calculate_skew(genome):
    # genome is a sequence of C,G,A,T
    skews = [0]
    for pair_value in genome:
        skew_value = -1 if pair_value == 'C' else (1 if pair_value == 'G' else 0) 
        skews.append(skews[-1] + skew_value)
    return skews

def minimum_skew(genome):
    # assume genome is from 5' -> 3' direction
    # the minimum is when we switch from decreasing to increasing, equals the ori
    skews = calculate_skew(genome)
    min_skew = min(skews)
    min_indices = [i for i, x in enumerate(skews) if x == min_skew]    
    return min_indices

def hamming_distance(genome1, genome2):
    if len(genome1) != len(genome2):
        raise ValueError('the two inputs must be same length')
    mismatches = 0
    for i in range(len(genome1)):
        if genome1[i] != genome2[i]:
            mismatches += 1
    return mismatches

def approximate_matches(pattern, genome, d):
    # Input: Strings Pattern and Text along with an integer d.
    # Output: All starting positions where Pattern appears as a substring of Text with at most d mismatches.
    matches = []
    for i in range(len(genome) - len(pattern) + 1):  # len(genome)
        sliced_genome = genome[i:i+len(pattern)]
        hd = hamming_distance(sliced_genome, pattern)  # O(len(pattern))
        if hd <= d:
            matches.append(i)
    return matches

def approximate_matches_count(pattern, genome, d):
    approx_matches = approximate_matches(pattern, genome, d)
    return len(approx_matches)

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


def frequent_words_with_mismatches_build_map(genome, k, d):
    # k length of the mer
    freqMap = dict()
    for i in range(len(genome) - k + 1):
        pattern = genome[i:i+k]
        neighborhood = neighbors(pattern, d)
        for neighbor in neighborhood:
            if neighbor not in freqMap:
                freqMap[neighbor] = 0
            freqMap[neighbor] += 1
    return freqMap

def frequent_words_with_mismatches(genome, k, d):
    freqMap = frequent_words_with_mismatches_build_map(genome, k, d)
    max_freq = max([freqMap[k] for k in freqMap])
    return [k for k in freqMap if freqMap[k] == max_freq]
    

def complement(genome):
    complementary_pairs = {
        'A' : 'T',
        'T' : 'A',
        'G' : 'C',
        'C' : 'G',
    }
    return ''.join([complementary_pairs[x] for x in genome][::-1])

def frequent_words_with_mismatches_plus_complement(genome, k, d):
    freqMap1 = frequent_words_with_mismatches_build_map(genome, k, d)
    freqMap2 = frequent_words_with_mismatches_build_map(complement(genome), k, d)

    mergedFreqMap = {k:freqMap1[k] for k in freqMap1}
    for k in freqMap2:
        if k in mergedFreqMap:
            mergedFreqMap[k] += freqMap2[k]
        else:
            mergedFreqMap[k] = freqMap2[k]
            

    max_freq = max([mergedFreqMap[k] for k in mergedFreqMap])
    return [k for k in mergedFreqMap if mergedFreqMap[k] == max_freq]