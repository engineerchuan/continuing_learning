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