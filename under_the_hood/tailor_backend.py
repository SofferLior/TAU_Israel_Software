def find_amino_seq(file_location, gene_code):
    amino_seq = ''
    file = open(file_location)
    found_it = False
    line = file.readline()
    while line != '':
        if line[0] == '>':
            if not found_it:
                gene_data = line.split('|')
                current_gene_code = gene_data[0][1:gene_data[0].find(' ')]
                if current_gene_code == gene_code:
                    found_it = True
            else:
                return amino_seq
        else:
            if found_it:
                amino_seq += line[:-1]
        line = file.readline()
    return amino_seq


def search_the_db(file_location, host_or_virus_name, search_by_virus, exact_match):
    file = open(file_location)
    line = file.readline()
    ans_lists = []
    while line != '':
        if line[0] == '>':
            gene_data = line.split('|')
            gene_code = gene_data[0][1:gene_data[0].find(' ')]
            virus_name = gene_data[0][gene_data[0].find(' ')+1:]
            gene_name = gene_data[1]
            host_name = gene_data[2]
            virus_lineage = gene_data[3]
            host_lineage = gene_data[4]
            if (search_by_virus) and (exact_match) and (virus_name == host_or_virus_name):
                ans_lists.append([gene_code,gene_name,virus_name,host_name,virus_lineage,host_lineage])
            elif (search_by_virus) and (not exact_match) and (host_or_virus_name in virus_name):
                ans_lists.append([gene_code,gene_name,virus_name,host_name,virus_lineage,host_lineage])
            elif (not search_by_virus) and (exact_match) and (host_or_virus_name == host_name):
                ans_lists.append([gene_code,gene_name,virus_name,host_name,virus_lineage,host_lineage])
            elif (not search_by_virus) and (not exact_match) and (host_or_virus_name in host_name):
                ans_lists.append([gene_code,gene_name,virus_name,host_name,virus_lineage,host_lineage])
        line = file.readline()
    return ans_lists