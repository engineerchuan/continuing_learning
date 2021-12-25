# encoding=utf8
# read in the data
import codecs
with codecs.open('dataset_2_6.txt', encoding='utf-8') as fid:    
    text = fid.readline().strip()
    pattern = fid.readline().strip()
    fid.close()
    
    count = 0
    for i in range(len(text) - len(pattern)+1):
        if text[i:i+len(pattern)] == pattern:
            count += 1
    print(count)