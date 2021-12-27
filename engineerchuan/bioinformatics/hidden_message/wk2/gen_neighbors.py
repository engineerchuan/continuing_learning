from lib_wk2 import neighbors

# encoding=utf8
# read in the data
import codecs
with codecs.open('data/dataset_3014_4.txt', encoding='utf-8') as fid:   
#if 1==1:
    #pattern = 'GAGG'
    genome = fid.readline().strip()#'ACG'
    d = int(fid.readline())
    print(' '.join(neighbors(genome, d)))