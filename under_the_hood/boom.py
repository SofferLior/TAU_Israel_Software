import csv
from Bio import SeqIO


def boom_normalized_words_score(seq, words_dict,word_size):
    score = 0
    for i in range(len(seq)-(word_size-1)):
        if seq[i:i+word_size] in words_dict:
            score = score + 1
    return score/(max(len(seq)-(word_size-1), len(words_dict)))


def all_genes_words_grades_table(file_location, genes_info_dict, genes_set, words_dict, write,word_size):
    if write:
        genes_table = open(file_location, mode='w', newline='')
        genes_table_writer = csv.writer(genes_table, delimiter=',', quoting=csv.QUOTE_MINIMAL)
        genes_table_writer.writerow(['grade', 'gene_code', 'virus_name', 'gene_name', 'host_name', 'virus_lineage', 'host_lineage'])
    grades_dict = {}
    for gene in genes_set:
        grade = boom_normalized_words_score(genes_info_dict[gene].protein[:200], words_dict.copy(), word_size)
        if write:
            genes_table_writer.writerow([grade, genes_info_dict[gene].gene_code, genes_info_dict[gene].virus_name, genes_info_dict[gene].gene_name, genes_info_dict[gene].host_name, genes_info_dict[gene].virus_lineage, genes_info_dict[gene].host_lineage])
        if grade not in grades_dict:
            grades_dict[grade] = set()
        grades_dict[grade].add(gene)
    return grades_dict


def get_words_dict_of_seq(seq, word_size):
    words_dict = set()
    for i in range(len(seq)-(word_size-1)):
        words_dict.add(str(seq[i:i+word_size]))
    return words_dict


def get_top_heuristic_grades(grades_dict, threshold):
    grades_list = list(grades_dict.keys())
    grades_list.sort(reverse=True)
    top_grades_set = set()
    i = 0
    while len(top_grades_set)< threshold:
        top_grades_set = top_grades_set.union(grades_dict[grades_list[i]])
        i = i + 1
    return top_grades_set


def get_boom_results(genes_info_dict, how_many_to_return):
    prf15 = SeqIO.read('C:/Users/ofeks/Documents/University/IGEM/bioinformatics/prf15.fasta', "fasta")[0:164]
    word_size = 3
    words_set = get_words_dict_of_seq(prf15.seq, word_size)
    grades_dict = all_genes_words_grades_table('C:/Users/ofeks/Documents/University/IGEM/bioinformatics/genes_words_table.csv', genes_info_dict,genes_info_dict.keys(), words_set, True, word_size)
    return get_top_heuristic_grades(grades_dict, how_many_to_return)
