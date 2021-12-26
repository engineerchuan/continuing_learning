from lib_wk2 import minimum_skew


# encoding=utf8
# read in the data
import codecs
with codecs.open('data/dataset_7_10.txt', encoding='utf-8') as fid:   
    genome = fid.readline().strip() 
    print(minimum_skew(genome))