from lib_wk2 import approximate_matches

# encoding=utf8
# read in the data
import codecs
with codecs.open('data/dataset_9_4.txt', encoding='utf-8') as fid:   
#if 1==1:
    #pattern = 'AAAAA'
    #genome = 'AACAAGCTGATAAACATTTAAAGAG'
    #d = 2
    pattern = fid.readline().strip() 
    genome = fid.readline().strip() 
    d = int(fid.readline().strip())
#genome1 = "GGGCCGTTGGT"
#genome2 = "GGACCGTTGAC"
    print(' '.join([str(x) for x in approximate_matches(pattern, genome, d)]))