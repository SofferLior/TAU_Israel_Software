from Bio import pairwise2
from Bio import SeqIO
from Bio.SubsMat.MatrixInfo import blosum62
from Bio.Seq import Seq


def similar_to_prf15_beginning(seq, ret_only_score, *parameters):
    r2_gene = SeqIO.read('C:/Users/ofeks/Documents/University/IGEM/bioinformatics/prf15.fasta', "fasta")
    if len(parameters) == 0:
        alignments = dp_beginning_alignment_with_extension(r2_gene.seq[:164], Seq(seq), blosum62, -11, -1, ret_only_score)
    elif len(parameters) == 3:
        alignments = dp_beginning_alignment_with_extension(r2_gene.seq[:164], Seq(seq), parameters[0], parameters[1], parameters[2], ret_only_score)
    else:
        print('wrong number of parameters, should get only 3')
    return alignments


def similar_to_prf16(seq, ret_only_score):
    r2_gene = SeqIO.read('C:/Users/ofeks/Documents/University/IGEM/bioinformatics/prf16.fasta', "fasta")
    alignments = pairwise2.align.globalds(r2_gene.seq, Seq(seq), blosum62, -11, -1, score_only=ret_only_score)
    return alignments


def dp_beginning_alignment_recursion(i, j, dp_dict, beginning, t, gap, similarity_matrix, beginning_length):
    if i == 0 and j == 0:
        return 0
    elif i == 0:
        return dp_dict[i, j-1] + gap
    elif j == 0:
        return dp_dict[i-1, j] + gap
    elif i == (beginning_length-1):
        val = max(dp_dict[i, j-1], dp_dict[i-1, j] + gap, dp_dict[i-1, j-1] + similarity_matrix[beginning[i], t[j]])
        if val > dp_dict['max_value']:
            dp_dict['max_value'] = val
            dp_dict['max_location'] = j
        return val
    else:
        return max(dp_dict[i, j-1] + gap, dp_dict[i-1, j] + gap, dp_dict[i-1, j-1] + similarity_matrix[beginning[i], t[j]])


def dp_beginning_alignment(beginning, other_seq, similarity_matrix, gap, score_only):
    dp_dict = {'max_value':float("-inf")}
    beginning = '_' + beginning
    other_seq = '_' + other_seq
    beginning_length = len(beginning)
    for aa in list(similarity_matrix.keys()):
        similarity_matrix[aa[1], aa[0]] = similarity_matrix[aa[0], aa[1]]
    for i in range(beginning_length):
        for j in range(len(other_seq)):
            dp_dict[i,j] = dp_beginning_alignment_recursion(i, j, dp_dict, beginning, other_seq, gap, similarity_matrix, beginning_length)
    if score_only:
        return dp_dict[len(beginning)-1, len(other_seq)-1]
    else:
        return pairwise2.align.globalds(beginning[1:], other_seq[1:dp_dict['max_location']+1], similarity_matrix, gap, gap, score_only=False)


def similarity_percentage(s, t):
    counter = 0
    for i in range(min(len(s),len(t))):
        if s[i] == t[i]:
            counter = counter + 1
    return counter/max(len(s),len(t))


def dp_beginning_alignment_with_extension(beginning, other_seq, similarity_matrix, gap, gap_extension, score_only):
    dp_dict = {'max_value':float("-inf")}
    beginning = '_' + beginning
    other_seq = '_' + other_seq
    beginning_length = len(beginning)
    for aa in list(similarity_matrix.keys()):
        similarity_matrix[aa[1], aa[0]] = similarity_matrix[aa[0], aa[1]]
    for i in range(beginning_length):
        for j in range(len(other_seq)):
            dp_dict[1, i, j] = dp_beginning_alignment_recursion_extension_diagonal(i, j, dp_dict, beginning, other_seq, similarity_matrix, beginning_length)
            dp_dict[2, i, j] = dp_beginning_alignment_recursion_extension_up(i, j, dp_dict, gap, gap_extension, beginning_length)
            dp_dict[3, i, j] = dp_beginning_alignment_recursion_extension_left(i, j, dp_dict, gap, gap_extension, beginning_length)
    if score_only:
        return max(dp_dict[1,len(beginning) - 1, len(other_seq) - 1], dp_dict[2,len(beginning) - 1, len(other_seq) - 1], dp_dict[3,len(beginning) - 1, len(other_seq) - 1])
    else:
        return pairwise2.align.globalds(beginning[1:], other_seq[1:dp_dict['max_location'] + 1], similarity_matrix, gap, gap_extension, score_only=False)


def dp_beginning_alignment_recursion_extension_left(i, j, dp_dict, gap, gap_extension, beginning_length):
    if j == 0:
        return float("-inf")
    elif i == (beginning_length-1):
        val = max(dp_dict[1,i, j-1], dp_dict[2,i, j-1], dp_dict[3,i, j-1])
        if val > dp_dict['max_value']:
            dp_dict['max_value'] = val
            dp_dict['max_location'] = j
        return val
    else:
        return max(dp_dict[1,i, j-1] + gap, dp_dict[2,i, j-1] + gap, dp_dict[3,i, j-1] + gap_extension)


def dp_beginning_alignment_recursion_extension_up(i, j, dp_dict, gap, gap_extension, beginning_length):
    if i == 0:
        return float("-inf")
    else:
        val = max(dp_dict[1, i-1, j] + gap, dp_dict[2, i-1, j] + gap_extension, dp_dict[3, i-1, j] + gap)
        if i == (beginning_length - 1) and val > dp_dict['max_value']:
            dp_dict['max_value'] = val
            dp_dict['max_location'] = j
        return val


def dp_beginning_alignment_recursion_extension_diagonal(i, j, dp_dict, beginning, other_seq, similarity_matrix, beginning_length):
    if i == 0 and j == 0:
        return 0
    elif i == 0:
        return float("-inf")
    elif j == 0:
        return float("-inf")
    else:
        val = similarity_matrix[beginning[i], other_seq[j]] + max(dp_dict[1, i-1, j - 1] , dp_dict[2, i-1, j - 1] , dp_dict[3, i-1, j - 1])
        if i == (beginning_length - 1) and val > dp_dict['max_value']:
            dp_dict['max_value'] = val
            dp_dict['max_location'] = j
        return val