import codecs, copy

def eulerian_cycle(adjacency_map):
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
            #node = adjacency_map.keys()[0]
            node = 3
            current_path = []
            # pick an outgoing edge
        else:
            # we finished one subcycle but need to continue with another
            for i in range(1, len(current_path) - 1):
                start_candidate = current_path[i]
                if len(unused_edges[start_candidate]) > 0:
                    node = start_candidate
                    current_path = current_path[i:] + current_path[1:i]
                    #unused_edges = copy.deepcopy(adjacency_map)
                    break
            
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

exercise_01_01_a()