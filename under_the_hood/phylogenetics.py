from Bio.Phylo.TreeConstruction import *
from Bio.Phylo import draw_ascii
from Bio import pairwise2, Phylo
from Bio.Seq import Seq
from io import StringIO
from ete3 import NCBITaxa, Tree
from Bio.SubsMat.MatrixInfo import *
from random import shuffle, randint
import math
from scipy.spatial import distance
import pickle
import os


def create_names_seqs_lists(genes_info_dict, genes_codes):
    names_list = []
    sequences_list = []
    for gene in genes_codes:
        names_list.append(genes_info_dict[gene].virus_name)
        sequences_list.append(genes_info_dict[gene].protein)
    return names_list, sequences_list


def random_sequance(seq):
    rand_seq = list(seq)
    shuffle(rand_seq)
    return ''.join(rand_seq)


def identical_seq_score(seq, scoring_matrix):
    score = 0
    for letter in seq:
        score = score + scoring_matrix[letter, letter]
    return score


def new_find_distance_between_seqs(seq1, seq2, scoring_matrix, gap, gap_extension):
    seq1_iden = identical_seq_score(seq1, scoring_matrix)
    seq2_iden = identical_seq_score(seq2, scoring_matrix)
    dist_between = pairwise2.align.globalds(Seq(seq1), Seq(seq2), scoring_matrix, gap, gap_extension, score_only=True)
    return (seq1_iden + seq2_iden - 2*dist_between)**0.5


def find_distance_between_seqs(seq1, seq2, scoring_matrix, gap, gap_extension): #using Feng at ell methode (1985)
    s_real = pairwise2.align.globalds(Seq(seq1), Seq(seq2), scoring_matrix, gap, gap_extension, score_only=True)
    rand_seq1 = random_sequance(seq1)
    rand_seq2 = random_sequance(seq2)
    s_rand = pairwise2.align.globalds(Seq(rand_seq1), Seq(rand_seq2), scoring_matrix, gap, gap_extension, score_only=True)
    s_iden = min(identical_seq_score(seq1, scoring_matrix), identical_seq_score(seq2, scoring_matrix))
    if s_real - s_rand > 0:
        s_eff = ((s_real - s_rand)/(s_iden - s_rand))
    else:
        s_eff = (0.001/(s_iden - s_rand)) # correction for the rare case that s_rand>s_iden as used in Da-Fei Feng, Russell F. Doolittle 1996
    return -(math.log(s_eff))*100


def create_distance_matrix(names_list, sequences_list, scoring_matrix, gap, gap_extension):
    distance_matrix = []
    for i in range(len(names_list)):
        distance_matrix.append([])
        for j in range(i):
            distance_matrix[i].append(new_find_distance_between_seqs(sequences_list[i], sequences_list[j], scoring_matrix, gap, gap_extension))
        distance_matrix[i].append(0)
    return DistanceMatrix(names_list, distance_matrix)


def create_tree_from_dm(dm):
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)
    return tree


def print_tree(tree):
    draw_ascii(tree)


def phylo_tree_2_ete_tree(tree,names_list):
    ncbi = NCBITaxa()
    handle = StringIO()
    Phylo.write(tree, handle, 'newick')
    newick_tree = handle.getvalue()
    name2taxid_dic = ncbi.get_name_translator(names_list)
    for i in range(len(names_list)):
        newick_tree = newick_tree.replace(names_list[i], str(name2taxid_dic[names_list[i]][0]))
    newick_tree = newick_tree.replace('\'', '')
    return Tree(newick_tree,format=1)


def get_ncbi_taxonomy_species_tree(names_list):
    ncbi = NCBITaxa()
    name2taxid_dic = ncbi.get_name_translator(names_list)
    taxid_list = []
    for i in range(len(names_list)):
        taxid_list.append(name2taxid_dic[names_list[i]][0])
    return ncbi.get_topology(taxid_list)


def compare_tree(tree1, tree2):
    return Tree.compare(tree1, tree2, unrooted=True)['rf']


def tree_builder(names_list, sequences_list, scoring_matrix, gap, gap_extension):
    dm = create_distance_matrix(names_list, sequences_list, scoring_matrix, gap, gap_extension)
    tree = create_tree_from_dm(dm)
    return phylo_tree_2_ete_tree(tree, names_list)


def find_closest_point(points_list, center_point):
    min_dist = float("inf")
    for point in points_list:
        dist = distance.euclidean(point, center_point)
        if dist < min_dist:
            min_dist = dist
            closest_point = point
    return closest_point


def hill_climbing_optimization(genes_info_dict, genes_codes):
    names_list, sequences_list = create_names_seqs_lists(genes_info_dict, genes_codes)
    ncbi_tree = get_ncbi_taxonomy_species_tree(names_list)
    matrixes = [blosum30, blosum35, blosum40, blosum45, blosum50, blosum55, blosum60, blosum62, blosum65, blosum70, blosum75, blosum80, blosum85, blosum90, blosum95, blosum100]
    memoization_climbing = {'min_result': float("inf")}
    #hill_climbing_recursion(memoization_climbing, ncbi_tree, names_list, sequences_list, matrixes, 7, -11, -1, float("inf"))
    hill_climbing_recursion_run_on_server(memoization_climbing, ncbi_tree, names_list, sequences_list, matrixes, 7, -11, -1, float("inf"))
    min_res = memoization_climbing['min_result']
    min_list = memoization_climbing['min_location']
    print(min_list)
    closest_point = find_closest_point(min_list, (7, -11, -1))
    return matrixes[closest_point[0]], closest_point[1], closest_point[2]


def hill_climbing_recursion(memoization_climbing, ncbi_tree, names_list, sequences_list, matrixes, matrix_num, original_gap, original_gap_extension, started_min):
    print("starting with center of matrix = {0}, gap = {1}, extension = {2} while min is {3}".format(matrix_num, original_gap, original_gap_extension, started_min))
    current_min = started_min
    for i in [matrix_num-1, matrix_num, matrix_num+1]:
        for gap_extension in [original_gap_extension-0.5, original_gap_extension, original_gap_extension+0.5]:
            for gap in [original_gap-1, original_gap, original_gap+1]:
                print(memoization_climbing)
                print('tring with blosum matrix {0} gap = {1} and extension = {2}'.format(i, gap, gap_extension))
                if (gap < 0) and (gap_extension < 0) and (i >= 0):
                    print("post inspections")
                    scoring_matrix = matrixes[i]
                    if (i, gap, gap_extension) not in memoization_climbing:
                        print('tree')
                        tree = tree_builder(names_list, sequences_list, scoring_matrix, gap, gap_extension)
                        print('dist')
                        dist = compare_tree(tree, ncbi_tree)
                        print('the result is {0}'.format(dist))
                        memoization_climbing[(i, gap, gap_extension)] = dist
                        print('saved')
                        if dist < memoization_climbing['min_result']:
                            print("updating global min result")
                            memoization_climbing['min_result'] = dist
                            memoization_climbing['min_location'] = [(i, gap, gap_extension)]
                        elif dist == memoization_climbing['min_result']:
                            memoization_climbing['min_location'].append((i, gap, gap_extension))
                        if dist < current_min:
                            print("updating local min result")
                            current_min = dist
                            continue_to = [(i, gap, gap_extension)]
                        elif (dist == current_min) and (dist < started_min):
                            continue_to.append((i, gap, gap_extension))
                    else:
                        print(memoization_climbing)
    print('got here')
    if started_min > current_min:
        started_min = current_min
        for parameters in continue_to:
            hill_climbing_recursion(memoization_climbing, ncbi_tree, names_list, sequences_list, matrixes, parameters[0], parameters[1], parameters[2], started_min)


def hill_climbing_recursion_run_on_server(memoization_climbing, ncbi_tree, names_list, sequences_list, matrixes, matrix_num, original_gap, original_gap_extension, started_min):
    current_min = started_min
    f1 = open('SERVER_DATA/store_ncbi_tree', 'wb')
    f2 = open('SERVER_DATA/store_names_list', 'wb')
    f3 = open('SERVER_DATA/store_sequences_list', 'wb')
    f4 = open('SERVER_DATA/store_matrixes', 'wb')
    pickle.dump(ncbi_tree, f1)
    pickle.dump(names_list, f2)
    pickle.dump(sequences_list, f3)
    pickle.dump(matrixes, f4)
    f1.close()
    f2.close()
    f3.close()
    f4.close()
    fscript = open('CODE/run_all_py_jobs.sh', 'w', newline='\n')
    fscript.write('#!/bin/bash\n')
    fscript.write('cd /tamir1/ofeks/igem/CODE/\n')
    for i in [matrix_num-1, matrix_num, matrix_num+1]:
        for gap_extension in [original_gap_extension-0.5, original_gap_extension, original_gap_extension+0.5]:
            for gap in [original_gap-1, original_gap, original_gap+1]:
                if (gap < 0) and (gap_extension < 0) and (i >= 0):
                    if (i, gap, gap_extension) not in memoization_climbing:
                        fscript.write('qsub -q tamirs2 -e ../JOBS/ERR/ -o ../JOBS/OUT/ -l cput=5:00:00,pmem=3gb,mem=3gb,pvmem=4gb,vmem=4gb ../JOBS/job_num{0},{1},{2}.sh\n'.format(i, gap, gap_extension))
                        f = open('JOBS/job_num{0},{1},{2}.sh'.format(i, gap, gap_extension), 'w', newline='\n')
                        f.write('#!/bin/bash\n')
                        f.write('cd /tamir1/ofeks/igem/CODE/\n')
                        f.write('module load python/python-anaconda3.6.5\n')
                        f.write('python run_on_server.py ''find_tree_and_dist_from_ncbi'' {0} {1} {2}'.format(i, gap, gap_extension))
                        f.close()
    fscript.close()
    input("Press Enter to continue...")
    continue_to = []
    for i in [matrix_num - 1, matrix_num, matrix_num + 1]:
        for gap_extension in [original_gap_extension - 0.5, original_gap_extension, original_gap_extension + 0.5]:
            for gap in [original_gap - 1, original_gap, original_gap + 1]:
                if (gap < 0) and (gap_extension < 0) and (i >= 0):
                    if (i, gap, gap_extension) not in memoization_climbing:
                        while not os.path.isfile('RESULTS/{0},{1},{2}'.format(i, float(gap), float(gap_extension))):
                            input("file RESULTS/{0},{1},{2} was not found, press enter to try again".format(i, float(gap), float(gap_extension)))
                        f = open('RESULTS/{0},{1},{2}'.format(i, float(gap), float(gap_extension)), 'r')
                        dist = int(float(f.readlines()[0]))
                        f.close()
                        if dist < memoization_climbing['min_result']:
                                print("updating global min result from {0} to {1}".format(memoization_climbing['min_result'], dist))
                                memoization_climbing['min_result'] = dist
                                memoization_climbing['min_location'] = [(i, gap, gap_extension)]
                        elif dist == memoization_climbing['min_result']:
                            memoization_climbing['min_location'].append((i, gap, gap_extension))
                        if dist < current_min:
                            print("updating local min result from {0} to {1}".format(current_min, dist))
                            current_min = dist
                            continue_to = [(i, gap, gap_extension)]
                        elif (dist == current_min) and (dist < started_min):
                            continue_to.append((i, gap, gap_extension))
    print('min is {0}, continuing to {1} locations'.format(memoization_climbing['min_result'], len(continue_to)))
    if started_min > current_min:
        started_min = current_min
        for parameters in continue_to:
            hill_climbing_recursion_run_on_server(memoization_climbing, ncbi_tree, names_list, sequences_list, matrixes, parameters[0], parameters[1], parameters[2], started_min)


