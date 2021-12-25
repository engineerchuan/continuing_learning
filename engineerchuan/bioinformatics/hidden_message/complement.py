# encoding=utf8
# read in the data
import codecs
complementary_pairs = {
    'A' : 'T',
    'T' : 'A',
    'G' : 'C',
    'C' : 'G',
}
with codecs.open('dataset_3_2.txt', encoding='utf-8') as fid:    
    input_line = fid.readline().strip()
    fid.close()

    print(''.join([complementary_pairs[x] for x in input_line][::-1]))

