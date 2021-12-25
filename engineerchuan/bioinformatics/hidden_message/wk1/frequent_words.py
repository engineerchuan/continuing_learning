# encoding=utf8
# read in the data
import codecs
with codecs.open('dataset_2_13.txt', encoding='utf-8') as fid:    
    puzzle = fid.readline().strip()
    n = int(fid.readline().strip())
    fid.close()

    def print_top_non_decoded_letters(n=1):
        char_freq_map = dict()
        for i in range(len(puzzle) - n + 1):
            candidate = puzzle[i:i+n]
            if candidate not in char_freq_map:
                char_freq_map[candidate] = 0
            char_freq_map[candidate] += 1

        # switch to tuples
        char_freq_tuples = [(letter, char_freq_map[letter]) for letter in char_freq_map]
        # filter by those not sorted
        max_freq = max([freq for (letter, freq) in char_freq_tuples])
        for (candidate, freq) in char_freq_tuples:
            if freq == max_freq:
                print('%10s %d' % (candidate, freq))

    print_top_non_decoded_letters(n=n)