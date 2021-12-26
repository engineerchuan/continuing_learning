from lib_wk2 import approximate_matches_count

# encoding=utf8
# read in the data
import codecs
with codecs.open('data/dataset_9_6.txt', encoding='utf-8') as fid:   
#if 1==1:
    #pattern = 'GAGG'
    #genome = 'TTTAGAGCCTTCAGAGG'
    #d = 2
    pattern = fid.readline().strip() 
    genome = fid.readline().strip() 
    d = int(fid.readline().strip())
#genome1 = "GGGCCGTTGGT"
#genome2 = "GGACCGTTGAC"
    print(approximate_matches_count(pattern, genome, d))