import csv
import scoring
from Bio.SubsMat.MatrixInfo import blosum62


class Gene:
    def __init__(self, gene_code, virus_name,gene_name,host_name,virus_lineage,host_lineage,sample_type,taxonomic_identifier,protein):
        self.gene_code = gene_code
        self.virus_name = virus_name
        self.gene_name = gene_name
        self.host_name = host_name
        self.virus_lineage = virus_lineage
        self.host_lineage = host_lineage
        self.sample_type = sample_type
        self.taxonomic_identifier = taxonomic_identifier
        self.protein = protein


def fasta_db_parser(file_location):
    file = open(file_location)
    line = file.readline()
    first = 1
    genes_info_dict = {}
    host_dict = {}
    virus_dict = {}
    tail_fiber_set = set()
    tail_spike_set = set()
    while line != '':
        if line[0] == '>':
            if first == 1:
                first = 0
            else:
                amino_seq = amino_seq.replace('\n', '')
                genes_info_dict[gene_code] = Gene(gene_code, virus_name,gene_name,host_name,virus_lineage,host_lineage,sample_type,taxonomic_identifier,amino_seq)
            gene_data = line.split('|')
            gene_code = gene_data[0][1:gene_data[0].find(' ')]
            virus_name = gene_data[0][gene_data[0].find(' ')+1:]
            gene_name = gene_data[1]
            host_name = gene_data[2]
            virus_lineage = gene_data[3]
            host_lineage = gene_data[4]
            sample_type = gene_data[5]
            taxonomic_identifier = gene_data[6]
            amino_seq = ''
            if 'tail fiber' in gene_name:
                tail_fiber_set.add(gene_code)
            if 'tail spike' in gene_name:
                tail_spike_set.add(gene_code)
            if host_name not in host_dict:
                host_dict[host_name] = set()
            host_dict[host_name].add(gene_code)
            if virus_name not in virus_dict:
                virus_dict[virus_name] = set()
            virus_dict[virus_name].add(gene_code)
        else:
            amino_seq = amino_seq + line
        line = file.readline()
    return genes_info_dict, virus_dict, host_dict, tail_fiber_set, tail_spike_set


def get_grades_of_genes(genes_info_dict, genes_set, similarity_matrix, gap, gap_extension, get_similarity):
    genes_grades = {}
    if get_similarity:
        genes_similarity = {}
    for gene in genes_set:
        if get_similarity:
            results = scoring.similar_to_prf15_beginning(genes_info_dict[gene].protein, False, similarity_matrix, gap, gap_extension)
            grade = results[0][2]
            similarity_percentage = scoring.similarity_percentage(results[0][0], results[0][1])
            genes_similarity[gene] = similarity_percentage
        else:
            grade = scoring.similar_to_prf15_beginning(genes_info_dict[gene].protein, True, similarity_matrix, gap, gap_extension)
        genes_grades[gene] = grade
    if get_similarity:
        return genes_grades, genes_similarity
    else:
        return genes_grades


def all_genes_grades_table_writer(file_location, genes_info_dict, genes_grades, genes_similarity):
    genes_table = open(file_location, mode='w', newline='')
    genes_table_writer = csv.writer(genes_table, delimiter=',', quoting=csv.QUOTE_MINIMAL)
    genes_table_writer.writerow(['grade', 'similarity_percentage','gene_code', 'virus_name', 'gene_name', 'host_name', 'virus_lineage', 'host_lineage'])
    for gene in genes_grades.keys():
        genes_table_writer.writerow([genes_grades[gene], genes_similarity[gene], genes_info_dict[gene].gene_code, genes_info_dict[gene].virus_name, genes_info_dict[gene].gene_name, genes_info_dict[gene].host_name, genes_info_dict[gene].virus_lineage, genes_info_dict[gene].host_lineage])


def only_highest_grade_per_virus(genes_info_dict, genes_grades):
    highest_gene = {}
    for gene in list(genes_grades.keys()):
        if genes_info_dict[gene].virus_name not in highest_gene.keys():
            highest_gene[genes_info_dict[gene].virus_name] = gene
        else:
            if genes_grades[highest_gene[genes_info_dict[gene].virus_name]] > genes_grades[gene]:
                genes_grades.pop(gene)
            else:
                genes_grades.pop(highest_gene[genes_info_dict[gene].virus_name])
                highest_gene[genes_info_dict[gene].virus_name] = gene


def only_top_grades(genes_grades, genes_info_dict, top_grades_threshold, differnet_taxonomy_threshold):
    grades = {}
    max_grade = float("-inf")
    min_grade = float("inf")
    sorted_out_genes_score = {}
    sorted_out_genes_similarity = {}
    in_list_taxonomy = set()
    for gene in genes_grades.keys(): # finding max score
        if max_grade < genes_grades[gene]:
            max_grade = genes_grades[gene]
        if min_grade < genes_grades[gene]:
            min_grade = genes_grades[gene]
        if genes_grades[gene] in grades:
            grades[genes_grades[gene]] = grades[genes_grades[gene]] + 1
        else:
            grades[genes_grades[gene]] = 1
    counter = 0
    min_threshold_grade = max_grade
    while counter < top_grades_threshold: # finding min score according to thrushold
        if min_threshold_grade in grades:
            counter = counter + grades[min_threshold_grade]
        min_threshold_grade = min_threshold_grade - 1
    for gene in list(genes_grades.keys()): # sorting only genes above the min score
        if genes_grades[gene] < min_threshold_grade:
            sorted_out_genes_score[gene] = genes_grades[gene]
            genes_grades.pop(gene)
        else:
            in_list_taxonomy.add(genes_info_dict[gene].virus_lineage)
    for gene in list(sorted_out_genes_score.keys()):
        if differnet_taxonomy_threshold == 0:
            return
        if genes_info_dict[gene].virus_lineage not in in_list_taxonomy:
            differnet_taxonomy_threshold = differnet_taxonomy_threshold - 1
            genes_grades[gene] = sorted_out_genes_score[gene]


def find_closest_to_prf16(file_location, genes_info_dict, genes_set):
    genes_table = open(file_location, mode='w', newline='')
    genes_table_writer = csv.writer(genes_table, delimiter=',', quoting=csv.QUOTE_MINIMAL)
    genes_table_writer.writerow(['grade', 'gene_code', 'virus_name', 'gene_name', 'host_name', 'virus_lineage', 'host_lineage'])
    max_similarity = float("-inf")
    for gene in genes_set:
        grade = scoring.similar_to_prf16(genes_info_dict[gene].protein, True)
        genes_table_writer.writerow([grade, genes_info_dict[gene].gene_code, genes_info_dict[gene].virus_name, genes_info_dict[gene].gene_name, genes_info_dict[gene].host_name, genes_info_dict[gene].virus_lineage, genes_info_dict[gene].host_lineage])
        if grade > max_similarity:
            max_similarity = grade
            max_gene = gene
    return max_gene
