from lib_wk2 import frequent_words_with_mismatches

# encoding=utf8
# read in the data
import codecs
with codecs.open('data/dataset_9_9.txt', encoding='utf-8') as fid:   
#if 1==1:
    #pattern = 'GAGG'
    genome = 'AACAAGCTGATAAACATTTAAAGAG'
    k = 5
    d = 1
    #pattern = fid.readline().strip() 
    genome = fid.readline().strip() 
    __ = fid.readline().strip().split()
    k = int(__[0])
    d = int(__[1])
    #d = int(fid.readline().strip())
#genome1 = "GGGCCGTTGGT"
#genome2 = "GGACCGTTGAC"
    print(frequent_words_with_mismatches(genome, k, d))