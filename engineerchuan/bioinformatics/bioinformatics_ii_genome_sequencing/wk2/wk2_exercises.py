import codecs, copy
import sys
sys.path.append("../wk1")
from wk1_exercises import graph_from_k_mers

def eulerian_cycle(adjacency_map, start = None):
    # input should be 
    # map[node1] = [node2, node3]
    current_path = None
    unused_edges = copy.deepcopy(adjacency_map)
    # check that it is balanced
    # check that it is strongly connected

    while True:
        #print("STARTING MEGA LOOP")
        # if everything is covered we are done!
        if sum([len(unused_edges[node]) for node in unused_edges.keys()]) == 0:
            break
        elif current_path is None:            
            if start is None:
                node = adjacency_map.keys()[0]
            else:
                node = start
            current_path = []
            # pick an outgoing edge
        else:
            # we finished one subcycle but need to continue with another
            print(current_path)
            for i in range(1, len(current_path) - 1):
                start_candidate = current_path[i]
                if len(unused_edges[start_candidate]) > 0:
                    node = start_candidate
                    current_path = current_path[i:] + current_path[1:i]
                    #unused_edges = copy.deepcopy(adjacency_map)
                    break
            print(node)
            
            
        # step 1, pick an outgoing edge
        while True:
            #print(current_path)
            #print(unused_edges)
            if len(unused_edges[node]) == 0:
                # we have reached a cycle
                current_path.append(node)
                break
            else:
                current_path.append(node)
                node = unused_edges[node].pop()       

    return current_path


def exercise_01_01_a():
    adjacency_map = {0:[3], 1:[0], 2:[1,6], 3:[2], 4:[2], 5:[4], 6:[5,8], 7:[9], 8:[7], 9:[6]}

    with codecs.open('data/dataset_203_2.txt', encoding='utf-8') as fid:

        adjacency_map = dict()
        line = fid.readline()
        while line != "":
            tuples = line.split(' -> ')
            key_ = int(tuples[0])
            values_ = [int(x) for x in tuples[1].split(',')]
            adjacency_map[key_] = values_
            line = fid.readline()   
    cycle = eulerian_cycle(adjacency_map)

    print('->'.join([str(x) for x in cycle]))


def eulerian_path(adjacency_map):
    # i will build upon eulerian_cycle
    # first compute the in degrees and outdegrees of each node.
    # find the two that are inbalanced, and draw the edge
    # then find the cycle starting with the inbalanced node.
    # then lop off the last one
    in_degree_out_degree = dict()
    # set the out degrees
    for key in adjacency_map:
        if key not in in_degree_out_degree:
            in_degree_out_degree[key] = [0, 0]
        in_degree_out_degree[key][1] += len(adjacency_map[key])
    # set the in degree
    for key in adjacency_map:
        values = adjacency_map[key]
        for value in values:
            if value not in in_degree_out_degree:
                in_degree_out_degree[value] = [0, 0]
            in_degree_out_degree[value][0] += 1
    # find the one node with in degree of 0
    #print(in_degree_out_degree)
    origin = [x for x in in_degree_out_degree if in_degree_out_degree[x][0] == (in_degree_out_degree[x][1]-1) ][0]
    #print(origin)
    # find the one node with out degree of 0
    terminus = [x for x in in_degree_out_degree if in_degree_out_degree[x][1] == (in_degree_out_degree[x][0] -1) ][0]
    # add the fake link
    if terminus not in adjacency_map:
        adjacency_map[terminus] = []
    adjacency_map[terminus].append(origin)
    cycle = eulerian_cycle(adjacency_map, start = None)
    # find the origin and then
    path = cycle[:-1]
    istart = path.index(origin)
    return path[istart:] + path[:istart]

    
def exercise_01_01_b():
    if True:
        adjacency_map = {0:[2], 1:[3], 2:[1], 3:[0,4], 6:[3,7], 7:[8], 8:[9], 9:[6]}

    with codecs.open('data/dataset_203_6.txt', encoding='utf-8') as fid:

        adjacency_map = dict()
        line = fid.readline()
        while line != "":
            tuples = line.split(' -> ')
            key_ = int(tuples[0])
            values_ = [int(x) for x in tuples[1].split(',')]
            adjacency_map[key_] = values_
            line = fid.readline()   
        cycle = eulerian_path(adjacency_map)

    print('->'.join([str(x) for x in cycle]))

def string_reconstruction(patterns):
    graph = graph_from_k_mers(patterns)
    path = eulerian_path(graph)
    return path[0] + ''.join([x[-1] for x in path[1:]])
    

def exercise_01_01_c():
    if True:
        patterns = ['CTTA','ACCA','TACC','GGCT','GCTT','TTAC']
    with codecs.open('data/dataset_203_7.txt', encoding='utf-8') as fid:
        patterns = []
        line = fid.readline()
        line = fid.readline()
        while line != "":
            patterns.append(line.strip())
            line = fid.readline()           
        print(string_reconstruction(patterns))

def solve_universal_string(k):
    import itertools
    all_possible = map(''.join, itertools.product('01', repeat=k))
    graph = graph_from_k_mers(all_possible)
    cycle = eulerian_cycle(graph)
    # adjust for the fact it is a circular string
    return cycle[0] + ''.join([x[-1] for x in cycle[1:-(k-1)]])

def generate_paired_k_mer_composition(dna, k, d):
    results = []
    for i in range(len(dna) - k - k - d + 1):
        first = dna[i:i+k]
        second = dna[i+k+d:i+k+d+k]
        results.append((first, second))
    results.sort()
    return results

#print(solve_universal_string(9))
#print(len('0000110010111101'))
#results = generate_paired_k_mer_composition('TAATGCCATGGGATGTT', 3, 2)
#print(' '.join(['(%s|%s)' % (first, second) for (first, second) in results]))

def graph_from_k_mer_pairs(kmers):
    # these are pairs separated by |
    k = (len(kmers[0]) - 1) // 2
    adjacency = dict()
    for sliced in kmers:
        first_k_minus_one_mer = sliced[:k-1] + '|' + sliced[k+1:2*k] 
        second_k_minus_one_mer = sliced[1:k] + '|' + sliced[k+2:] 

        if first_k_minus_one_mer not in adjacency:
            adjacency[first_k_minus_one_mer] = []
        adjacency[first_k_minus_one_mer].append(second_k_minus_one_mer)
    return adjacency


def string_reconstruction_k_d(patterns, k, d):
    graph = graph_from_k_mer_pairs(patterns)
    #print(graph)
    path = eulerian_path(graph)
    k = (len(path[0]) - 1) // 2
    firsts = [x[:k] for x in path]
    seconds = [x[k+1:] for x in path]
    firsts_unrolled = firsts[0] + ''.join([x[-1] for x in firsts[1:]])
    seconds_unrolled = seconds[0] + ''.join([x[-1] for x in seconds[1:]])
    #print(firsts_unrolled)
    #print(seconds_unrolled)
    return firsts_unrolled[:k+d+1] + seconds_unrolled


def exercise_01_01_d():
    patterns = ['GAGA|TTGA','TCGT|GATG','CGTG|ATGT','TGGT|TGAG','GTGA|TGTT','GTGG|GTGA','TGAG|GTTG','GGTC|GAGA','GTCG|AGAT']

    k = 4
    d = 2
    with codecs.open('data/dataset_204_16.txt', encoding='utf-8') as fid:
        patterns = []
        __ = fid.readline().split(' ')
        k = int(__[0])
        d = int(__[1])
        line = fid.readline()
        while line != "":
            patterns.append(line.strip())
            line = fid.readline()           
        print(string_reconstruction_k_d(patterns, k, d))
    
exercise_01_01_d()
    