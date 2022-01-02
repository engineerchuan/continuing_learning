import codecs

codon_to_k = {
    'AAA':'K',
    'AAC':'N',
    'AAG':'K',
    'AAU':'N',
    'ACA':'T',
    'ACC':'T',
    'ACG':'T',
    'ACU':'T',
    'AGA':'R',
    'AGC':'S',
    'AGG':'R',
    'AGU':'S',
    'AUA':'I',
    'AUC':'I',
    'AUG':'M',
    'AUU':'I',
    'CAA':'Q',
    'CAC':'H',
    'CAG':'Q',
    'CAU':'H',
    'CCA':'P',
    'CCC':'P',
    'CCG':'P',
    'CCU':'P',
    'CGA':'R',
    'CGC':'R',
    'CGG':'R',
    'CGU':'R',
    'CUA':'L',
    'CUC':'L',
    'CUG':'L',
    'CUU':'L',
    'GAA':'E',
    'GAC':'D',
    'GAG':'E',
    'GAU':'D',
    'GCA':'A',
    'GCC':'A',
    'GCG':'A',
    'GCU':'A',
    'GGA':'G',
    'GGC':'G',
    'GGG':'G',
    'GGU':'G',
    'GUA':'V',
    'GUC':'V',
    'GUG':'V',
    'GUU':'V',
    'UAA':'*',
    'UAC':'Y',
    'UAG':'*',
    'UAU':'Y',
    'UCA':'S',
    'UCC':'S',
    'UCG':'S',
    'UCU':'S',
    'UGA':'*',
    'UGC':'C',
    'UGG':'W',
    'UGU':'C',
    'UUA':'L',
    'UUC':'F',
    'UUG':'L',
    'UUU':'F'
}

def translate_rna_to_amino_acid(rna_string):
    assert len(rna_string) % 3 == 0
    amino_acid_string = ''
    for i in range(0, len(rna_string), 3):
        codon = codon_to_k[rna_string[i:i+3]]
        if codon == '*':
            return amino_acid_string
        else:
            amino_acid_string += codon

def reverse_complement(dna):
    # dna is a string
    mapped = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    return ''.join([mapped[x] for x in dna][::-1])

def exercise_03_01_a():

    with codecs.open('data/dataset_96_4.txt', encoding='utf-8') as fid:

        rna_string = fid.readline().strip()
        print(translate_rna_to_amino_acid(rna_string))

def convert_rna_to_dna(rna):
    return rna.replace('U', 'T')

# This is way too slow
def encoding_dna_substrings(dna, amino_acids, rna_to_amino_map):
    # step 1, reverse the map
    amino_to_dna_map = dict()
    for rna_, amino_ in rna_to_amino_map.items():
      
        if amino_ not in amino_to_dna_map:
            amino_to_dna_map[amino_] = []
        amino_to_dna_map[amino_].append(convert_rna_to_dna(rna_))
    import itertools
    candidates = itertools.product(*[amino_to_dna_map[x] for x in amino_acids])
    #print("created map, candidates %d", len([x for x in candidates]))
    answers = []
    i = 0
    for candidate in candidates:
        i+=1
        if i % 100 == 0:
            print(i)
        joined = ''.join(candidate)
        if (joined in dna):
            answers.append(joined)
        if reverse_complement(joined) in dna:
            answers.append(reverse_complement(joined))
    return answers

def exercise_03_01_b():
    dna = 'ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA'
    amino_acids = 'MA'
    rna_to_amino_map = codon_to_k
    with codecs.open('data/dataset_96_7.txt', encoding='utf-8') as fid:

        dna = fid.readline().strip()
        amino_acids = fid.readline().strip()

        answers = encoding_dna_substrings(dna, amino_acids, rna_to_amino_map)
        print('\n'.join(answers))

#exercise_03_01_b()

#def encoding_dna_substrings_with_position(dna, amino_acids, rna_to_amino_map):
#    answers = encoding_dna_substrings(dna, amino_acids, rna_to_amino_map)
#   start_positions = set()
#    for answer in answers:

def exercise_03_01_c():
    amino_acids = 'VKLFPWFNGY'
    rna_to_amino_map = codon_to_k
    with codecs.open('data/Bacillus_brevis.txt', encoding='utf-8') as fid:
        dna = ''.join([x.strip() for x in fid.readlines()])
        answers = encoding_dna_substrings(dna, amino_acids, rna_to_amino_map)
        print('\n'.join(answers))

exercise_03_01_c()
