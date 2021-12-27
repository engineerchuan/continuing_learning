def longest_common_substring(string1, string2):
    # initiate the matrix, first index will be for string1
    # second index for string 2
    # each matrix element will be the longest length
    # and the index in string1
    mat = [[(0, "") for i in range(len(string2)+1)] for j in range(len(string1)+1)]
    for i in range(len(string1)):
        for j in range(len(string2)):
            if string1[i] == string2[j]:
                (max_diagonal, string_so_far) = mat[i][j]
                mat[i+1][j+1] = (max_diagonal+1, string_so_far + string1[i])
            else:
                mat[i+1][j+1] = (0, "")
    # find the max
    max_so_far = 0
    lcs_so_far = ""
    for i in range(len(mat)):
        for j in range(len(mat[0])):
            (max_len, lcs) = mat[i][j]
            if max_len > max_so_far:
                max_so_far = max_len
                lcs_so_far = lcs
    return (max_so_far, lcs_so_far)

def exercise_01_02_a():
    a1 = 'atgaccgggatactgataaaaaaaagggggggggcgtacacattagataaacgtatgaagtacgttagactcggcgccgccg'
    a2 = 'acccctattttttgagcagatttagtgacctggaaaaaaaatttgagtacaaaacttttccgaataaaaaaaaaggggggga'
    a3 = 'tgagtatccctgggatgacttaaaaaaaagggggggtgctctcccgatttttgaatatgtaggatcattcgccagggtccga'
    (len_1_2, lcs) = longest_common_substring(a1, a2)
    print(longest_common_substring(lcs, a3))


exercise_01_02_a()

#print(longest_common_substring("bcde", "bce"))
