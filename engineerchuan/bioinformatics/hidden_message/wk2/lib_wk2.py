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
    for i in range(len(genome) - len(pattern) + 1):
        sliced_genome = genome[i:i+len(pattern)]
        hd = hamming_distance(sliced_genome, pattern)
        if hd <= d:
            matches.append(i)
    return matches
