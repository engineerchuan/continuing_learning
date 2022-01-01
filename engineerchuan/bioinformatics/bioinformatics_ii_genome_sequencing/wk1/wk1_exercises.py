import codecs

def generate_k_mer_composition(dna, k):
    # generate the k mer composition
    # A, C, G, T
    composition = []
    for i in range(len(dna) - k + 1):
        composition.append(dna[i:i+k])
    # sort lexicographically
    composition.sort()
    return composition

def exercise_01_01_a():
    # for sample problem
    #if True:
    #    k = 5
    #    dna = 'CAATCCAAC'

    with codecs.open('data/dataset_197_3.txt', encoding='utf-8') as fid:
        k = int(fid.readline().strip())
        dna = fid.readline().strip()
        print('\n'.join(generate_k_mer_composition(dna, k)))


def string_from_genome_path(dnas):
    # input is list of strings
    n_len = len(dnas[0])
    string = [dnas[0]]
    for i in range(1, len(dnas)):
        assert len(dnas[i]) == len(dnas[i-1])
        assert dnas[i][:(n_len-1)] == dnas[i-1][1:n_len]
        string.append(dnas[i][-1])
    return ''.join(string)

def exercise_01_01_b():
    # for sample problem
    #if True:
    #    dnas = ['ACCGA', 'CCGAA', 'CGAAG', 'GAAGC', 'AAGCT']

    with codecs.open('data/dataset_198_3.txt', encoding='utf-8') as fid:
        dnas = []
        line = fid.readline()
        while line != "":
            dnas.append(line.strip())
            line = fid.readline()
        print('read in X lines')
        print(string_from_genome_path(dnas))


def overlap_graph_problem(k_mer_patterns):
    # input k mer patterns (all of same length)
    # output adjacency graph
    # handle formatting separately
    k = len(k_mer_patterns[0])
    prefix_to_kmer = dict()
    suffix_to_kmer = dict()
    for pattern in k_mer_patterns:
        # add to prefixes
        prefix = pattern[:k-1]
        if prefix not in prefix_to_kmer:
            prefix_to_kmer[prefix] = []
        prefix_to_kmer[prefix].append(pattern)

        # add to suffixes
        suffix = pattern[1:k]
        if suffix not in suffix_to_kmer:
            suffix_to_kmer[suffix] = []
        suffix_to_kmer[suffix].append(pattern)

    # create the adjacency graph
    adjacency_graph = dict()
    for suffix in suffix_to_kmer:
        for pattern1 in suffix_to_kmer[suffix]:
            if suffix in prefix_to_kmer:
                if pattern1 not in adjacency_graph:
                    adjacency_graph[pattern1] = []
                adjacency_graph[pattern1] += prefix_to_kmer[suffix]
    return adjacency_graph

def pretty_print_overlap_graph(k_mer_patterns):
    adjacency_graph = overlap_graph_problem(k_mer_patterns)
    for node1 in adjacency_graph:
        print("%s -> %s" % (node1, ",".join(adjacency_graph[node1])))


def exercise_01_01_c():
    #if True:
    #    k_mer_patterns = ['ATGCG', 'GCATG', 'CATGC', 'AGGCA', 'GGCAT', 'GGCAC']
        #overlap_graph_problem(k_mer_patterns)
    with codecs.open('data/dataset_198_10.txt', encoding='utf-8') as fid:
        k_mer_patterns = []
        line = fid.readline()
        while line != "":
            k_mer_patterns.append(line.strip())
            line = fid.readline()
        pretty_print_overlap_graph(k_mer_patterns)


def construct_universal_k_string(k):
    # use brute force method
    # 1.) Generate all possible strings
    import itertools
    all_possible = map(''.join, itertools.product('01', repeat=k))

    # breadth first search
    def helper(constructed, left):
        #print(constructed, left)
        if len(left) == 0:
            return [constructed]
        elif constructed == '':
            subanswers = []
            for candidate in left:
                __ = helper(candidate, [x for x in left if x != candidate])
                subanswers += __
            return subanswers
        else:
            subanswers = []
            for candidate in left:
                #print("*"*80)
                if constructed[(-k+1):] == candidate[:(k-1)]:
                    #print("recursive path")
                    __ = helper(constructed + candidate[-1], [x for x in left if x != candidate])
                    subanswers += [x for x in __ if len(x) == (len(all_possible) + (k-1))]
            return subanswers
    return helper('', all_possible)

def path_graph(text, k):
    # return adjacency list
    adjacency = dict()
    for i in range(len(text) - k + 1):
        sliced = text[i:i+k]
        first_k_minus_one_mer = sliced[:-1]
        second_k_minus_one_mer = sliced[1:]

        if first_k_minus_one_mer not in adjacency:
            adjacency[first_k_minus_one_mer] = []
        adjacency[first_k_minus_one_mer].append(second_k_minus_one_mer)
    return adjacency

def pretty_print_path_graph(text, k):
    adjacency_graph = path_graph(text, k)
    for node1 in adjacency_graph:
        print("%s -> %s" % (node1, ",".join(adjacency_graph[node1])))


def exercise_01_01_d():
    if True:
        k = 4
        text = 'AAGATTCTCTAAGA'
        #overlap_graph_problem(k_mer_patterns)
    with codecs.open('data/dataset_199_6.txt', encoding='utf-8') as fid:
        k = int(fid.readline().strip())
        text = fid.readline().strip()
       
        pretty_print_path_graph(text, k)


def exercise_01_01_e():
    # question 1
    # no you could not reconstruct it
    # question 2
    print("2 mer")
    pretty_print_path_graph('TAATGCCATGGGATGTT', 2)
    print("3 mer")
    pretty_print_path_graph('TAATGCCATGGGATGTT', 3)
    print("4 mer")
    pretty_print_path_graph('TAATGCCATGGGATGTT', 4)
    # question 3
    pretty_print_path_graph('TAATGCCATGGGATGTT', 3)
    pretty_print_path_graph('TAATGGGATGCCATGTT', 3)
    # they are the same


def graph_from_k_mers(kmers):
    # a collection of patterns of length k
    # a de brujin graph adjacency list
    k = len(kmers[0])
    #adjacency = overlap_graph_problem(kmers)
    #print(adjacency)
    #de_brujin_graph = dict()
    #for first in adjacency:
    #    if first[:(k-1)] not in de_brujin_graph:
    ##        de_brujin_graph[first[:(k-1)]] = []
    #    de_brujin_graph[first[:(k-1)]].append(first[1:])#

    #for second in adjacency[first]:
    #        if second[:(k-1)] not in de_brujin_graph:
    #            de_brujin_graph[second[:(k-1)]] = []
    #        de_brujin_graph[second[:(k-1)]].append(second[1:])                
    #for k in de_brujin_graph:
    #   de_brujin_graph[k] = list(set(de_brujin_graph[k]))
    #return de_brujin_graph
    adjacency = dict()
    for sliced in kmers:
        first_k_minus_one_mer = sliced[:-1]
        second_k_minus_one_mer = sliced[1:]

        if first_k_minus_one_mer not in adjacency:
            adjacency[first_k_minus_one_mer] = []
        adjacency[first_k_minus_one_mer].append(second_k_minus_one_mer)
    return adjacency

def pretty_print_graph_from_k_mers(kmers):
    adjacency_graph = graph_from_k_mers(kmers)
    for node1 in adjacency_graph:
        print("%s -> %s" % (node1, ",".join(adjacency_graph[node1])))    

def exercise_01_01_f():
    #if True:
    #    kmers = ['GAGG','CAGG','GGGG','GGGA','CAGG','AGGG','GGAG']
    with codecs.open('data/dataset_200_8.txt', encoding='utf-8') as fid:
        kmers = []
        line = fid.readline()
        while line != "":
            kmers.append(line.strip())
            line = fid.readline()        
      
        pretty_print_graph_from_k_mers(kmers)

print(len(construct_universal_k_string(3)))