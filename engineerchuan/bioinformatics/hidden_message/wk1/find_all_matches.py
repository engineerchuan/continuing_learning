# encoding=utf8
# read in the data
import codecs
with codecs.open('../Vibrio_cholerae.txt', encoding='utf-8') as fid:    
    #pattern = fid.readline().strip()
    genome = fid.readline().strip()
    pattern = 'CTTGATCAT'
    fid.close()
    answers = []
    for i_start in range(len(genome) - len(pattern) + 1):
        candidate = genome[i_start:i_start+ len(pattern)]
        #print(candidate)
        #print(pattern)
        if candidate == pattern:
            answers.append(i_start)
    print(' '.join([str(x) for x in answers]))

