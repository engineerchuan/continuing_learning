# encoding=utf8
# read in the data
import codecs
with codecs.open('E_coli.txt') as fid:    
    genome = fid.readline().strip()
    print('read in genome of length %d' % len(genome))
    #params = fid.readline().strip().split(' ')
    #k = int(params[0])
    #l = int(params[1])
    #t = int(params[2])
    k = 9
    t = 3
    l = 500
    char_freq_map = dict()
    
    def print_top_non_decoded_letters(genome, char_freq_map, i_slice_start, l=50, k=1, t=0, build_map_everytime=False):
        answers = []
        n = k
        input_slice = genome[i_slice_start:i_slice_start+l]

        if build_map_everytime or i_slice_start == 0:
            # build it up naively
            for i in range(len(input_slice) - n + 1):
                candidate = input_slice[i:i+n]
                if candidate not in char_freq_map:
                    char_freq_map[candidate] = 0
                char_freq_map[candidate] += 1
        else:
            last_candidate = input_slice[i_slice_start-1:i_slice_start-1+n]
            print(char_freq_map)
            char_freq_map[last_candidate] -=1
            next_candidate = input_slice[i_slice_start+l-1:i_slice_start+l-1+n]
            char_freq_map[next_candidate] +=1
            #raise Exception("not implemented")

        # switch to tuples
        char_freq_tuples = [(letter, char_freq_map[letter]) for letter in char_freq_map]
        # filter by those not sorted
        #max_freq = max([freq for (letter, freq) in char_freq_tuples])
        for (candidate, freq) in char_freq_tuples:
            if freq >= t:
                answers.append(candidate)
                #print('%10s %d' % (candidate, freq))
        return answers
    
    total_answers = []
    for i_start in range(len(genome) - l + 1):
        #print(i_start)
        if i_start % 10000 == 0:
            print(i_start)
        input_genome_slice = genome[i_start:i_start+l]    
        answers = print_top_non_decoded_letters(genome = genome, char_freq_map=char_freq_map, i_slice_start = i_start, l=l, k=k, t=t)
        total_answers += answers
    print(' '.join(list(set(total_answers))))