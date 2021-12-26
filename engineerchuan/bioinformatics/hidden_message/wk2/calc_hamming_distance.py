from lib_wk2 import hamming_distance

# encoding=utf8
# read in the data
import codecs
with codecs.open('data/dataset_9_3.txt', encoding='utf-8') as fid:   
    genome1 = fid.readline().strip() 
    genome2 = fid.readline().strip()
#genome1 = "GGGCCGTTGGT"
#genome2 = "GGACCGTTGAC"
    print(hamming_distance(genome1, genome2))