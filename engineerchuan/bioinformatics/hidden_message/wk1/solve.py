# encoding=utf8
# read in the data
import codecs
with codecs.open('hidden_message.txt', encoding='utf-8') as fid:    
    puzzle = fid.read().replace('\n', '')
    fid.close()

    decoded_mapping = dict()

    def replace_letter(puzzle, coded_letter, decoded_letter):
        if coded_letter in decoded_mapping:
            raise Exception("letter already in the mapping")
        decoded_mapping[coded_letter] = decoded_letter
        return puzzle.replace(coded_letter, decoded_letter)

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
        char_freq_tuples = [(letter, letter_freq) for (letter, letter_freq) in char_freq_tuples if letter not in decoded_mapping]
        char_freq_tuples.sort(lambda x,y: y[1] - x[1])
        print("*"*80)
        print("top 10 strings of length %d" % n)
        print("*"*80)
        for (candidate, freq) in char_freq_tuples[:10]:            
            print('%10s %d' % (candidate, freq))
            

    
    puzzle = replace_letter(puzzle, u'8', 'E')
    puzzle = replace_letter(puzzle, u'4', 'H')
    puzzle = replace_letter(puzzle, u';', 'T')
    #puzzle = replace_letter(puzzle, u'5', 'N')
    #puzzle = replace_letter(puzzle, u'‡', 'L')
    puzzle = replace_letter(puzzle, u'†', 'R')
    #puzzle = replace_letter(puzzle, u')', 'S')
    #replace_letter(u'\u2021', 'I')
    #replace_letter(u')', 'O')
    print_top_non_decoded_letters(n=1)
    print_top_non_decoded_letters(n=2)
    print_top_non_decoded_letters(n=3)
    print_top_non_decoded_letters(n=4)
    print(puzzle)