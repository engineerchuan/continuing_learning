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

def generate_theoretical_spectrum(amino_acid_string, integer_mass_table, linear=False):
    # integer_mass_table is a map
    # amino_acid_string
    n = len(amino_acid_string)
    subpeptides = [[]]
    if not linear:
        for length in range(1, n):
            for i in range(0, n):
                if (i + length) > n:
                    subpeptides.append(amino_acid_string[i:] +amino_acid_string[:((i+length) % n)])
                else:
                    subpeptides.append(amino_acid_string[i:(i+length)])
    else:
        for length in range(1, n):
            for i in range(0, n):
                if (i + length) < n:
                    subpeptides.append(amino_acid_string[i:(i+length)])

    subpeptides.append(amino_acid_string)
    masses = []
    for subpeptide in subpeptides:
        masses.append(sum([integer_mass_table[x] for x in subpeptide]))
    masses.sort()
    return masses

def generate_cyclospectrum(masses):
    # masses are all the individual masses
    n = len(masses)
    subpeptides = [0]
    doubler = masses + masses
    for length in range(1, n):
        for i in range(0, n):
            subpeptides.append(sum(doubler[i:i+length]))
    subpeptides.append(sum(masses))
    subpeptides.sort()
    return subpeptides


def generate_linearspectrum(masses):
    # masses are all the individual masses
    n = len(masses)
    subpeptides = [0]
    for length in range(1, n):
        for i in range(0, n):
            if i + length < n:
                subpeptides.append(sum(masses[i:i+length]))
    subpeptides.append(sum(masses))
    subpeptides.sort()
    return subpeptides


def exercise_03_02_d():
    integer_mass_table = dict()

    with codecs.open('data/integer_mass_table.txt', encoding='utf-8') as fid:

        line = fid.readline().strip()
        while line != "":
            __ = line.strip().split()
            integer_mass_table[__[0]] = int(__[1])
            line = fid.readline().strip()
        fid.close()
    with codecs.open('data/dataset_98_4.txt', encoding='utf-8') as fid:
        peptide = fid.readline().strip()

        masses = generate_theoretical_spectrum(peptide, integer_mass_table)
        print(' '.join([str(x) for x in masses]))
    

def branch_and_bound(integer_mass_table, masses):
    # convert integer mass table into a list
    unique_masses = list(set(integer_mass_table.values()))
    #candidate_peptides = [tuple()]
    candidate_peptides = [()]
    parent_mass = masses[-1]
    final_peptides = []
    # take all the candidate peptides
    while len(candidate_peptides) > 0:
        #print(candidate_peptides)
        new_candidate_peptides = []
        for candidate in candidate_peptides:
            for mass in unique_masses:
                new_candidate_peptides.append(candidate + (mass,))
        candidate_peptides = new_candidate_peptides
        #print(candidate_peptides)
        survivors = []
        for candidate in candidate_peptides:
            candidate_spectrum = generate_cyclospectrum(candidate)
            linear_spectrum = generate_linearspectrum(candidate)
            
            all_in = min([x in masses for x in linear_spectrum])
            #print(candidate, linear_spectrum, all_in)
            if sum(candidate) == parent_mass:
                if list(candidate_spectrum) == masses:
                    final_peptides.append(candidate)
            elif sum(candidate) > parent_mass:
                pass
            elif all_in == False:
                pass
            else:
                survivors.append(candidate)
            #print("len survivors", len(survivors))
        candidate_peptides = survivors
    return final_peptides

import sys
def exercise_03_02_e():
    integer_mass_table = dict()
    with codecs.open('data/integer_mass_table.txt', encoding='utf-8') as fid:

        line = fid.readline().strip()
        while line != "":
            __ = line.strip().split()
            integer_mass_table[__[0]] = int(__[1])
            line = fid.readline().strip()
        fid.close()   

    with codecs.open('data/dataset_100_6.txt', encoding='utf-8') as fid:
        masses = [int(x) for x in fid.readline().strip().split(' ')]
        print(masses)

        #masses = generate_theoretical_spectrum(peptide, integer_mass_table)

        #masses = [0, 113, 128, 186, 241, 299, 314, 427]
        #masses = [0, 113, 128, 186]
        final_peptides = branch_and_bound(integer_mass_table, masses)
        for final in final_peptides:
            sys.stdout.write('-'.join([str(x) for x in final]) + ' ')
        print("")


            

def quiz():

    integer_mass_table = dict()
    with codecs.open('data/integer_mass_table.txt', encoding='utf-8') as fid:

        line = fid.readline().strip()
        while line != "":
            __ = line.strip().split()
            integer_mass_table[__[0]] = int(__[1])
            line = fid.readline().strip()
        fid.close()   

    desired = [0, 71, 101, 113, 131, 184, 202, 214, 232, 285, 303, 315, 345, 416]
    for candidate in ['TMLA', 'ALTM', 'IAMT', 'MTAI', 'TAIM', 'MIAT']:
        print(candidate, generate_theoretical_spectrum(candidate, integer_mass_table))
        if desired == generate_theoretical_spectrum(candidate, integer_mass_table):
            print(candidate)

    desired = [0, 71, 99, 101, 103, 128, 129, 199, 200, 204, 227, 230, 231, 298, 303, 328, 330, 332, 333]
    for candidate in ['VAQ', 'CET', 'ETC', 'AQV', 'CTV', 'TCE']:
        print(candidate, generate_theoretical_spectrum(candidate, integer_mass_table, linear=True))
        if desired == generate_theoretical_spectrum(candidate, integer_mass_table, linear=True):
            print(candidate)

    dnas = ['CCUCGUACAGAAAUCAAC', 'CCAAGAACAGAUAUCAAU', 'CCAAGUACAGAGAUUAAC', 'CCGAGGACCGAAAUCAAC']
    for dna in dnas:
        amino_acids = 'PRTEIN'
        rna_to_amino_map = codon_to_k
        print(dna)
        answers = encoding_dna_substrings(dna, amino_acids, rna_to_amino_map)
        print(answers)


quiz()