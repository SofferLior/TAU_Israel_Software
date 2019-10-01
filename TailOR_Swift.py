import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
import csv
from PIL import Image, ImageTk
import ttkwidgets


###########3main functions###################33
from operator import itemgetter

class Gene:
    def __init__(self, gene_code, virus_name,gene_name,host_name,alignment_val,virus_lineage,host_lineage):
        self.gene_code = gene_code
        self.virus_name = virus_name
        self.gene_name = gene_name
        self.host_name = host_name
        self.alignment= alignment_val
        self.virus_lineage = virus_lineage
        self.host_lineage = host_lineage

def load_csv():
    file_location='top_grades+and_tails_alignment.csv'
    genes_file = open(file_location, mode='r', newline='')
    genes_file_reader = csv.reader(genes_file, delimiter=',')
    genes_info_dict = {}
    virus_dict={}
    host_dict={}
    next( genes_file_reader)
    for row in genes_file_reader:
        gene_code= row[2]
        virus_name=row[3]
        gene_name =row[4]
        host_name=row[5]
        virus_lineage=row[6]
        host_lineage=row[7]
        alignment_val=float(row[0])
        genes_info_dict[gene_code]= Gene(gene_code, virus_name,gene_name,host_name,alignment_val,virus_lineage,host_lineage)
        if virus_name not in virus_dict:
            virus_dict[virus_name] = set()
        virus_dict[virus_name].add(gene_code)
        if host_name not in host_dict:
            host_dict[host_name] = set()
        host_dict[host_name].add(gene_code)

    return genes_info_dict,virus_dict,host_dict


def perfect_match_finder(name,searchby,genes_info_dict,virus_dict,host_dict):
    result=[]
    if searchby=='Virus':
        if name in virus_dict.keys():
            for gene in virus_dict[name]:
                if float(genes_info_dict[gene].alignment)>0:
                    result.append([genes_info_dict[gene].alignment,genes_info_dict[gene].gene_code,genes_info_dict[gene].gene_name,genes_info_dict[gene].virus_name,genes_info_dict[gene].host_name,genes_info_dict[gene].virus_lineage,genes_info_dict[gene].host_lineage])

    else:
        if name in host_dict.keys():
            for gene in host_dict[name]:
                if float(genes_info_dict[gene].alignment)>0:
                    result.append([genes_info_dict[gene].alignment,genes_info_dict[gene].gene_code,genes_info_dict[gene].gene_name,genes_info_dict[gene].virus_name,genes_info_dict[gene].host_name,genes_info_dict[gene].virus_lineage,genes_info_dict[gene].host_lineage])

    result=sorted(result, key=itemgetter(0),reverse=True)
    return result


def search_the_db(host_or_virus_name, search_by_virus, exact_match):
    file_location = 'virushostdb.formatted.cds.faa'
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


def search_the_db_for_tails(host_or_virus_name, search_by_virus, exact_match):
    file_location = 'virushostdb.formatted.cds.faa'
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
            if (search_by_virus) and (exact_match) and (virus_name == host_or_virus_name) and ('tail' in gene_name):
                ans_lists.append([gene_code,gene_name,virus_name,host_name,virus_lineage,host_lineage])
            elif (search_by_virus) and (not exact_match) and (host_or_virus_name in virus_name) and ('tail' in gene_name):
                ans_lists.append([gene_code,gene_name,virus_name,host_name,virus_lineage,host_lineage])
            elif (not search_by_virus) and (exact_match) and (host_or_virus_name == host_name) and ('tail' in gene_name):
                ans_lists.append([gene_code,gene_name,virus_name,host_name,virus_lineage,host_lineage])
            elif (not search_by_virus) and (not exact_match) and (host_or_virus_name in host_name) and ('tail' in gene_name):
                ans_lists.append([gene_code,gene_name,virus_name,host_name,virus_lineage,host_lineage])
        line = file.readline()
    return ans_lists


def find_amino_seq(gene_code):
    file_location='virushostdb.formatted.cds.faa'
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


################## Fustion Functions -FindTargetSeq########################3
pyocin_seq = "MTTNTPKYGGLLTDIGAAALATASAAGKKWQPTHMLIGDAGGAPGDTPDPLPSAAQKSLINQRHRAQLNRLFVSDKNANTLVAEVVLPVEVGGFWIREIG" \
             "LQDADGKFVAVSNCPPSYKAAMESGSARTQTIRVNIALSGLENVQLLIDNGIIYATQDWVKEKVAADFKGRKILAGNGLLGGGDLSADRSIGLAPSGVTA" \
             "GSYRSVTVNANGVVTQGSNPTTLAGYAIGDAYTKADTDGKLAQKANKATTLAGYGITDALRVDGNAVSSSRLAAPRSLAASGDASWSVTFDGSANVSAPL" \
             "SLSATGVAAGSYPKVTVDTKGRVTAGMALAATDIPGLDASKLVSGVLAEQRLPVFARGLATAVSNSSDPNTATVPLMLTNHANGPVAGRYFYIQSMFYPD" \
             "QNGNASQIATSYNATSEMYVRVSYAANPSIREWLPWQRCDIGGSFTKEADGELPGGVNLDSMVTSGWWSQSFTAQAASGANYPIVRAGLLHVYAASSNFI" \
             "YQTYQAYDGESFYFRCRHSNTWFPWRRMWHGGDFNPSDYLLKSGFYWNALPGKPATFPPSAHNHDVGQLTSGILPLARGGVGSNTAAGARSTIGAGVPAT" \
             "ASLGASGWWRDNDTGLIRQWGQVTCPADADASITFPIPFPTLCLGGYANQTSAFHPGTDASTGFRGATTTTAVIRNGYFAQAVLSWEAFGR"
pyocin_len = 691


def is_domain_file(file_path):
    with open(file_path, 'r') as my_file:
        data = my_file.read()
        if data.startswith("hmmscan"):
            return True
        return False


def is_structure_file(file_path):
    with open(file_path, 'r') as my_file:
        data = my_file.read()
        if data.startswith("{\"feature\""):
            return True
        return False


def find_target_seq_from_file(target_domain_path):
    with open(target_domain_path, 'r') as my_file:
        data = my_file.read()
        return find_target_seq(data)


def find_target_seq(domain_file):
    seq_found = False
    seq = ""
    for line in domain_file.splitlines():
        if line.startswith("Query"):
            seq_found = True
        elif seq_found and line != "":
            seq += line
        elif seq_found and line == "":
            break
    seq = seq.replace(" ", "")
    return seq, len(seq)


def fusion_protein(target_domain_path, pyocin_ind, target_ind):
    target_seq, target_len = find_target_seq_from_file(target_domain_path)
    if pyocin_ind < 0 or pyocin_ind > pyocin_len - 1 or target_ind < 0 or target_ind > target_len - 1:
        return 0
    return pyocin_seq[:pyocin_ind] + target_seq[target_ind:]


############### SecondaryStructureParse ################################

pyocin_cuts = [(168, 204), (209, 483), (486, 488), (494, 498), (506, 510), (516, 526), (528, 545), (548, 607), (610, 615), (625, 631), (636, 642), (649, 661), (666, 670), (675, 682)]
pyocin_domains = [(4, 158, "DUF3751")]


def parse_secondary_from_file(file_path):
    with open(file_path, 'r') as my_file:
        data = my_file.read()
        return parse_secondary_from_seq(data)


def parse_secondary_from_seq(seq):
    data = seq.replace('\n', '').split("\"featureString\":\"")[1].split("\",\"reliability\"")[0].replace("\\n", "")
    res = []
    current_pos = 0
    curr_len = 0
    total_len = 0
    for pos in data:
        total_len += 1
        if pos != 'L' and curr_len == 0:
            current_pos += 1
        elif pos != 'L' and curr_len != 0:
            res.append((current_pos, current_pos + curr_len))
            current_pos += (curr_len + 1)
            curr_len = 0
        else:
            curr_len += 1
    return res, total_len


def parse_domain_from_file(file_path):
    with open(file_path, 'r') as my_file:
        data = my_file.read()
        return parse_domain_from_string(data)


def parse_domain_from_string(data):
    start_delim = "=" * 130
    domain_start = False
    first_blank = False
    res = []
    for line in data.splitlines():
        if line == start_delim:
            domain_start = True
        elif domain_start and line == "" and not first_blank:
            first_blank = True
        elif domain_start and first_blank and line == "":
            break
        elif domain_start and first_blank and line != "":
            domain_start = False
            first_blank = False
            arr = line.split(" ")
            res.append((int(arr[3]), int(arr[4]), arr[0]))
        else:
            continue
    return res


def cut_locations(domain_path, secondary_path):
    target_domain = parse_domain_from_file(domain_path)
    secondary, protein_len = parse_secondary_from_file(secondary_path)
    res = []
    for second in secondary:
        for ind, domain in enumerate(target_domain):
            if second[1] <= domain[0]:
                res.append(second)
                break
            elif (second[0] < domain[0]) and (second[1] >= domain[0]) and (second[1] < domain[1]):
                res.append((second[0], domain[0] - 1))
                break
            elif second[0] >= domain[1] and (ind + 1 == len(target_domain)):
                res.append(second)
                break
            elif (second[0] < domain[0]) and (second[1] > domain[1]):
                res.append((second[0], domain[0] - 1))
                res.append((domain[1], second[1]))
                break
    return res, protein_len, target_domain


def fusion(domain_path, secondary_path):
    target_cuts, protein_len, target_domain = cut_locations(domain_path, secondary_path)
    return target_cuts, protein_len, target_domain




###################


class StColors(object):
    orange       = '#f5be2e'
    bright_green = '#b7f731'
    dark_grey    = '#191919'
    mid_grey     = '#323232'
    light_grey   = '#c8c8c8'
    purple       = '#ab72cc'
    white        = '#FFFFFF'
    torkiz       = '#00CED1'
    dark_torkiz  = '#12a9ab'
    dark_white   = '#eae7ec'
    cyten        = "#e5f1f4"


def on_tab_selected(event):
    selected_tab = event.widget.select()
    tab_text = event.widget.tab(selected_tab, "text")


def perfect_match_finder_command():
    EmptyRow = ttk.Label(tab1)
    EmptyRow.grid(row=9, column=0, padx=15, pady=0)
    genes_info_dict, virus_dict, host_dict = load_csv()
    searchby=SearchBYEntry.get()
    name=NameEntry.get()
    if searchby=='' or name=='':
        return
    global ans_list
    ans_list=perfect_match_finder(name,searchby,genes_info_dict,virus_dict,host_dict)
    if len(ans_list)>0:

        create_result_in_gui_with_alignment(ans_list)
        options_for_continue(11,tab1)
    else:
        ErrorLabel2 = ttk.Label(tab1, text='We do not have data on the request')
        ErrorLabel2.grid(row=10, column=0, columnspan=2)
        NameEntry.delete(0, 'end')



def options_for_continue(next_row,tab):
    EmptyRow = ttk.Label(tab)
    EmptyRow.grid(row=next_row, column=0, padx=15, pady=0)
    NewSearchButton=ttk.Button(tab, text='New Search',command=new_search)
    ExportButton=ttk.Button(tab,text='Export Results',command=excel_export)
    NewSearchButton.grid(row=next_row+1, column=1)
    ExportButton.grid(row=next_row+1, column=2)

    EmptyRow = ttk.Label(tab)
    EmptyRow.grid(row=next_row+2, column=0, padx=15, pady=0)


def excel_export():
    if tab_parent.index('current') == 2:
        fasta_export_seq()
    else:
        excel_export_list()


def excel_export_list():
    searchby=SearchBYEntry.get()
    name=NameEntry.get()
    dash='-'
    filename=dash.join(('TailOR_Swift',searchby,name,'Finder'))
    location_for_save = filedialog.askdirectory()
    backslash = '/'
    file=backslash.join((location_for_save,filename))
    file=file+'.csv'
    genes_table = open(file, mode='w', newline='')
    genes_table_writer = csv.writer(genes_table, delimiter=',', quoting=csv.QUOTE_MINIMAL)
    if len(ans_list[0])==7:
        genes_table_writer.writerow(['alignment','gene_code', 'gene_name','virus_name', 'host_name', 'virus_lineage', 'host_lineage'])
    else:
        genes_table_writer.writerow(['gene_code', 'gene_name', 'virus_name', 'host_name', 'virus_lineage', 'host_lineage'])
    for gene in ans_list:
        genes_table_writer.writerow(gene)


def fasta_export_seq():
    if tab_parent.index('current') == 2:
        gene_code = GeneCodeEntry.get()
        seq_to_fasta=aa_seq
    else:
        gene_code="fusion_leg"
        seq_to_fasta=fusion_seq
    location_for_save=filedialog.askdirectory()
    file=location_for_save+'/'+ gene_code+'.fasta'
    txtfile = open(file, mode='w')
    firstrow='>'+ gene_code+'\n'
    txtfile.write(firstrow)
    for i in range(int(len(seq_to_fasta)/80)):
        txtfile.write(seq_to_fasta[i*80:((i+1)*80-1)])
    txtfile.close()


def new_search():
    tabidx=tab_parent.index('current')
    if tabidx==1:
        for x in tab1.grid_slaves():
            if int(x.grid_info()["row"])>8:
                x.grid_forget()
        SearchBYEntry.delete(0,'end')
        NameEntry.delete(0,'end')
    if tabidx==2:
        for x in tab2.grid_slaves():
            if int(x.grid_info()["row"])>8:
                x.grid_forget()
        GeneCodeEntry.delete(0, 'end')
    if tabidx==3:
        for x in tab3.grid_slaves():
            if int(x.grid_info()["row"])>1:
                x.grid_forget()
        create_fusion()

def create_result_in_gui_with_alignment(ans_list):
    for x in tab1.grid_slaves():
        if int(x.grid_info()["row"]) > 9:
            x.grid_forget()
    result_table = ttk.Treeview(tab1, style="mystyle.Treeview")

    result_table["columns"] = ("one", "two", "three", "four", "five","six")
    result_table.heading("#0", text="Alignment")
    result_table.heading("one", text="Gene Code")
    result_table.heading("two", text="Gene Name")
    result_table.heading("three", text="Virus Name")
    result_table.heading("four", text="Host Name")
    result_table.heading("five", text="Virus Lineage")
    result_table.heading("six", text="Host Lineage")

    result_table.column("#0", width=70)
    result_table.column("one", width=120)
    result_table.column("two", width=120)
    result_table.column("three", width=120)
    result_table.column("four", width=120)
    result_table.column("five", width=120)
    result_table.column("six", width=120)
    for ans in ans_list:
        result_table.insert("", "end", text=ans[0], values=(ans[1], ans[2], ans[3], ans[4], ans[5],ans[6]))

    ysb = ttk.Scrollbar(tab1, orient='vertical', command=result_table.yview)
    # xsb = ttk.Scrollbar(tab1, orient='horizontal', command=result_table.xview)
    result_table.grid(row=10, column=0, columnspan=3)
    ysb.grid(row=10, column=3, sticky='ns')
    # xsb.grid(row=11,column=0,columnspan=4,sticky='ew')
    result_table.configure(yscroll=ysb.set)


def create_result_in_gui_no_alignment(ans_list):
    for x in tab1.grid_slaves():
        if int(x.grid_info()["row"]) > 9:
            x.grid_forget()
    result_table = ttk.Treeview(tab1, style="mystyle.Treeview")

    result_table["columns"] = ("one", "two", "three", "four", "five")
    result_table.heading("#0", text="Gene Code")
    result_table.heading("one", text="Gene Name")
    result_table.heading("two", text="Virus Name")
    result_table.heading("three", text="Host Name")
    result_table.heading("four", text="Virus Lineage")
    result_table.heading("five", text="Host Lineage")

    result_table.column("#0", width=120)
    result_table.column("one", width=120)
    result_table.column("two", width=120)
    result_table.column("three", width=120)
    result_table.column("four", width=120)
    result_table.column("five", width=120)
    for ans in ans_list:
        result_table.insert("", "end", text=ans[0], values=(ans[1], ans[2], ans[3], ans[4], ans[5]))

    ysb = ttk.Scrollbar(tab1, orient='vertical', command=result_table.yview)
    # xsb = ttk.Scrollbar(tab1, orient='horizontal', command=result_table.xview)
    result_table.grid(row=10, column=0, columnspan=3)
    ysb.grid(row=10, column=3, sticky='ns')
    # xsb.grid(row=11,column=0,columnspan=4,sticky='ew')
    result_table.configure(yscroll=ysb.set)


def db_searcher():
    EmptyRow = ttk.Label(tab1)
    EmptyRow.grid(row=9, column=0, padx=15, pady=0)
    searchby=SearchBYEntry.get()
    name= NameEntry.get()
    if searchby=='' or name=='':
        return
    global ans_list
    ans_list=search_the_db(name,searchby=='Virus',False)
    if len(ans_list)>0:
        create_result_in_gui_no_alignment(ans_list)
        options_for_continue(11,tab1)
    else:
        ErrorLabel2 = ttk.Label(tab1, text='We do not have data on the request')
        ErrorLabel2.grid(row=10, column=0, columnspan=2)
        NameEntry.delete(0, 'end')


def tail_searcher():
    EmptyRow = ttk.Label(tab1)
    EmptyRow.grid(row=9, column=0, padx=15, pady=0)
    searchby=SearchBYEntry.get()
    name= NameEntry.get()
    if searchby=='' or name=='':
        return
    global ans_list
    ans_list=search_the_db_for_tails(name,searchby=='Virus',False)
    if len(ans_list)>0:
        create_result_in_gui_no_alignment(ans_list)
        options_for_continue(11,tab1)
    else:
        ErrorLabel2 = ttk.Label(tab1, text='We do not have data on the request')
        ErrorLabel2.grid(row=10, column=0, columnspan=2)
        NameEntry.delete(0, 'end')


def find_seq():
    global aa_seq
    gene_code=GeneCodeEntry.get()
    aa_seq=find_amino_seq(gene_code)
    if aa_seq=='':
        ErrorLabel2 = ttk.Label(tab2, text='We do not have data on the request')
        ErrorLabel2.grid(row=10, column=0, columnspan=2)
        GeneCodeEntry.delete(0, 'end')
    else:
        AAtext = tk.Text(tab2)
        AAtext.insert('1.0',aa_seq,('text'))
        AAtext.grid(row=10,column=0,padx=15, pady=10,columnspan=10)
        AAtext.tag_configure('text',font=("Calibri", 12), relief='raised')
        options_for_continue(11,tab2)


def create_fusion():
    global domains_filename
    domains_filename = tk.filedialog.askopenfilename(initialdir="/", title="Select Domains file",
                                            filetypes=(("txt files", "*.txt"), ("all files", "*.*")))
    secondary_structure_filename = tk.filedialog.askopenfilename(initialdir="/", title="Select Secondary Structure file",
                                                                 filetypes=
                                                                 (("txt files", "*.txt"), ("all files", "*.*")))

    if secondary_structure_filename!="" and domains_filename!="":
        if is_domain_file(domains_filename) and is_structure_file(secondary_structure_filename):
            for x in tab3.grid_slaves():
                if int(x.grid_info()["row"]) > 1:
                    x.grid_forget()

            global target_cuts
            global protein_len
            target_cuts, protein_len, target_domain=fusion(domains_filename,secondary_structure_filename)
            Heading=ttk.Label(tab3,text="Target Domains",font=('Calibri', '12','bold'))
            Heading.grid(row=4, column=5, padx=15, pady=10)
            Domain = ttk.Label(tab3, text="Name", font=('Calibri', '10','bold'), foreground=StColors.dark_grey)
            Domain.grid(row=4, column=7, columnspan=2,padx=15, pady=10)
            Start = ttk.Label(tab3, text="Start", font=('Calibri', '10','bold'), foreground=StColors.dark_grey)
            Start.grid(row=4, column=10, columnspan=2,padx=15, pady=10)
            End = ttk.Label(tab3, text="End", font=('Calibri', '10','bold'), foreground=StColors.dark_grey)
            End.grid(row=4, column=12, columnspan=2,padx=15, pady=10)

            for i in range(len(target_domain)):
                Domain = ttk.Label(tab3, text=target_domain[i][2], font=('Calibri', '10'),foreground=StColors.dark_grey)
                Domain.grid(row=5+i,column=7,columnspan=2)
                Start = ttk.Label(tab3, text=target_domain[i][0], font=('Calibri', '10'),foreground=StColors.dark_grey)
                Start.grid(row=5+i,column=10, columnspan=2)
                End = ttk.Label(tab3, text=target_domain[i][1], font=('Calibri', '10'),foreground=StColors.dark_grey)
                End.grid(row=5+i,column=12, columnspan=2)

            global nextrow
            nextrow=5+i

            EmptyRow=ttk.Label(tab3)
            EmptyRow.grid(row=nextrow,column=0, padx=15, pady=0)

            Cutting_Headling=ttk.Label(tab3,text="Available Cutting Location",font=('Calibri', '12','bold'))
            Cutting_Headling.grid(row=nextrow+1, column=5, columnspan=2,padx=15, pady=10)

            Pyocin_Domains_Headling=ttk.Label(tab3,text="Pyocin:",font=('Calibri', '12','bold'))
            Pyocin_Domains_Headling.grid(row=nextrow+2, column=5, padx=15, pady=10)

            global Pyocin_Domains_Table
            global Target_Domains_Table
            Pyocin_Domains_Table=ttkwidgets.CheckboxTreeview(tab3)
            for cut in pyocin_cuts:
                 Pyocin_Domains_Table.insert("", "end", text=str(cut[0])+'-'+str(cut[1]))
            Pyocin_Domains_Table.grid(row=nextrow+3, column=5)

            Target_Domains_Headling=ttk.Label(tab3,text="Target:",font=('Calibri', '12','bold'))
            Target_Domains_Headling.grid(row=nextrow+2, column=8, padx=15, pady=10)

            Target_Domains_Table=ttkwidgets.CheckboxTreeview(tab3)
            for cut in target_cuts:
                Target_Domains_Table.insert("", "end", text=str(cut[0])+'-'+str(cut[1]))
            Target_Domains_Table.grid(row=nextrow + 3, column=8)

            ysb = ttk.Scrollbar(tab3, orient='vertical', command=Pyocin_Domains_Table.yview)
            ysb.grid(row=nextrow+3, column=6, sticky='ns')
            Pyocin_Domains_Table.configure(yscroll=ysb.set)

            ysb = ttk.Scrollbar(tab3, orient='vertical', command=Target_Domains_Table.yview)
            ysb.grid(row=nextrow+3, column=9, sticky='ns')
            Target_Domains_Table.configure(yscroll=ysb.set)


            EmptyRow = ttk.Label(tab3)
            EmptyRow.grid(row=nextrow + 4, column=0, padx=15, pady=0)
            Continue_Label = ttk.Label(tab3,text="Please choose cutting location for the Pyocin and the target")
            Continue_Label.grid(row=nextrow + 5, column=5,columnspan=4, padx=15, pady=0)

            NewSearchButton = ttk.Button(tab3, text='New Search', command=new_search)
            ExportButton = ttk.Button(tab3, text='Get the Sequence', command=create_fusion_sequence)

            EmptyRow = ttk.Label(tab3)
            EmptyRow.grid(row=nextrow + 6, column=0, padx=15, pady=0)

            NewSearchButton.grid(row=0, column=12)
            ExportButton.grid(row=nextrow + 5, column=11)
        else:
            EmptyRow = ttk.Label(tab3)
            EmptyRow.grid(row=5, column=0, padx=15, pady=0)

            ErrorFiles = ttk.Label(tab3, text='Please Load Domains file first and Secondary Structure file second')
            ErrorFiles.grid(row=7, column=5, columnspan=5)

    else:
        EmptyRow=ttk.Label(tab3)
        EmptyRow.grid(row=5,column=0, padx=15, pady=0)

        ErrorFiles = ttk.Label(tab3, text='Please Load all files!')
        ErrorFiles.grid(row=7, column=5, columnspan=5)


def create_fusion_sequence():
    cutting_location_target=Target_Domains_Table.get_checked()
    cutting_location_pyocin=Pyocin_Domains_Table.get_checked()

    for x in tab3.grid_slaves():
        if int(x.grid_info()["row"]) >nextrow+6 :
            x.grid_forget()

    if len(cutting_location_target)>1 or len(cutting_location_pyocin)>1:
        EmptyRow=ttk.Label(tab3)
        EmptyRow.grid(row=nextrow+7,column=5, padx=15, pady=0)
        ErrorLocation = ttk.Label(tab3, text='Please select only one possible location')
        ErrorLocation.grid(row=nextrow+8, column=5, columnspan=5)
    elif len(cutting_location_target)==0 or len(cutting_location_pyocin)==0:
        EmptyRow=ttk.Label(tab3)
        EmptyRow.grid(row=nextrow+7,column=5, padx=15, pady=0)
        ErrorLocation = ttk.Label(tab3, text='Please select possible location')
        ErrorLocation.grid(row=nextrow+8, column=5, columnspan=5)
    else:
        for x in tab3.grid_slaves():
            if int(x.grid_info()["row"]) > nextrow + 4:
                x.grid_forget()

        target_tuple_cut=int(cutting_location_target[0][1:])-1
        from_cutting_location_target=target_cuts[target_tuple_cut][0]
        to_cutting_location_target=target_cuts[target_tuple_cut][1]

        pyocin_tuple_cut=int(cutting_location_pyocin[0][1:])-1
        from_cutting_location_pyocin=pyocin_cuts[pyocin_tuple_cut][0]
        to_cutting_location_pyocin=pyocin_cuts[pyocin_tuple_cut][1]

        global scaleentryTarget
        global scaleentryPyocin
        ScaleTarget=ttk.Label(tab3,text="Select specific location for Target")
        ScaleTarget.grid(row=nextrow+6,column=5, columnspan=2,padx=15, pady=0)
        scaleentryTarget = ttkwidgets.ScaleEntry(tab3, scalewidth=200, entrywidth=3, from_=from_cutting_location_target, to=to_cutting_location_target)
        scaleentryTarget.grid(row=nextrow+6,column=8, columnspan=2)

        ScalePyocin=ttk.Label(tab3,text="Select specific location for Pyocin")
        ScalePyocin.grid(row=nextrow+5,column=5, columnspan=2,padx=15, pady=0)
        scaleentryPyocin = ttkwidgets.ScaleEntry(tab3, scalewidth=200, entrywidth=3, from_=from_cutting_location_pyocin, to=to_cutting_location_pyocin)
        scaleentryPyocin.grid(row=nextrow+5,column=8, columnspan=2)

        EmptyRow=ttk.Label(tab3)
        EmptyRow.grid(row=nextrow+7,column=5, padx=15, pady=0)
        ContinueButton = ttk.Button(tab3, text='Continue', command=specific_location_selection)
        ContinueButton.grid(row=nextrow+8,column=5)

def specific_location_selection():
    pyocin_ind = scaleentryPyocin.value
    target_ind = scaleentryTarget.value
    global fusion_seq
    fusion_seq = fusion_protein(domains_filename, pyocin_ind, target_ind)
    for x in tab3.grid_slaves():
        if int(x.grid_info()["row"]) > nextrow + 4:
            x.grid_forget()
    if fusion_seq==0:
        ErrorLocationSelection = ttk.Label(tab3, text='Please enter valid values')
        ErrorLocationSelection.grid(row=nextrow+7, column=5, columnspan=5)
    else:
        for x in tab3.grid_slaves():
            if int(x.grid_info()["row"]) > 1:
                x.grid_forget()
        FusionLabel = ttk.Label(tab3, text='The sequence for your new requested leg is ready!')
        FusionLabel.grid(row=2, column=1, columnspan=10,padx=15)
        FusionLocationLabel = ttk.Label(tab3, text='Pyocin: 0 to '+str(pyocin_ind)+'   Target: '+str(target_ind)+' to '
                                                   + str(protein_len))
        FusionLocationLabel.grid(row=3, column=1, columnspan=10, padx=15)

        new_fusion_len=protein_len-target_ind+pyocin_ind
        FusionLenLabel = ttk.Label(tab3, text='Total length: '+str(new_fusion_len))
        FusionLenLabel.grid(row=4, column=1, columnspan=10, padx=15)

        AAtext = tk.Text(tab3)
        AAtext.insert('1.0',fusion_seq,('text'))
        AAtext.grid(row=5,column=0,padx=15, pady=10,columnspan=15)
        AAtext.tag_configure('text',font=("Calibri", 12), relief='raised')

        ExportButton = ttk.Button(tab3, text='Export Results', command=fasta_export_seq)
        ExportButton.grid(row=7,column=5)


def set_app_style():
    style = ttk.Style()
    style.theme_create("st_app", parent="alt", settings={
        ".": {"configure": {"background": StColors.white,
                            "relief": "flat",
                            "highlightcolor": StColors.dark_grey}},

        "TLabel": {"configure": {"foreground": StColors.purple,
                                 "padding": 0,
                                 "font": ("Calibri", 12)}},

        "TNotebook": {"configure": {"padding": 0,"background": StColors.torkiz }},
        "TNotebook.Tab": {"configure": {"padding": [20, 20],
                                        "foreground": "white", "background": StColors.torkiz, "bordercolor": StColors.torkiz, "font" : ('Calibri', '12', 'bold')},
                          "map": {"background": [("selected", StColors.dark_torkiz)],
                                  "bordercolor": [("selected", StColors.dark_torkiz)],
                                  "expand": [("selected", [1, 1, 1, 0])]}},

        "TCombobox": {"configure": {"selectbackground": "white",
                                    "fieldbackground": "white",
                                    "selectforeground": "black",
                                    "arrowcolor": "white",
                                    "arrowsize ": 16,
                                    "bordercolor ": "white",
                                    "background": StColors.mid_grey,
                                    "foreground": "black"}},

        "TButton": {"configure": {"font": ("Calibri", 12),
                                  "background": StColors.cyten,
                                  "foreground": "black", "justify": "center"},
                    "map": {"background": [("active", StColors.purple)],
                            "foreground": [("active", 'white')]}},

        "TEntry": {"configure": {"foreground": "black", "bordercolor ": "white"}},
        "Horizontal.TProgressbar": {"configure": {"background": StColors.mid_grey}},
        "TTreeView": {"configure":{"font": ("Calibri", 12) }}
    })
    style.theme_use("st_app")


gui = tk.Tk()
gui.title("TailOR Swift")
#gui.geometry("500x280")

set_app_style()

myvar = tk.IntVar()

style = ttk.Style(gui)
style.configure('lefttab.TNotebook', tabposition='wn')
style.configure("Tab", focuscolor=StColors.dark_torkiz)
style.configure('Custom.TLabel', foreground="black", font=('Calibri', '8'))


style.configure("mystyle.Treeview", highlightthickness=0, bd=0, font=('Calibri', 11)) # Modify the font of the body
style.configure("mystyle.Treeview.Heading", font=('Calibri', 12),foreground=StColors.purple) # Modify the font of the headings
style.layout("mystyle.Treeview", [('mystyle.Treeview.treearea', {'sticky': 'nswe'})]) # Remove the borders

#s=ttk.Style()
#s.theme_use('st_app')


tab_parent = ttk.Notebook(gui, style='lefttab.TNotebook',width=850, height=650)

#s.configure('.', font=('Tahoma', 12))
#s.configure("TNotebook", padding=6, relief="flat",background="#ccc")


tab0 = ttk.Frame(tab_parent, width=450, height=450)
tab1 = ttk.Frame(tab_parent, width=450, height=450)
tab2 = ttk.Frame(tab_parent, width=450, height=450)
tab3 = ttk.Frame(tab_parent, width=450, height=450)


tab_parent.bind("<<NotebookTabChanged>>", on_tab_selected)
tab_parent.add(tab0, text="          Home         ")
tab_parent.add(tab1, text="   Tail Finder        ")
tab_parent.add(tab2, text="Get AA sequence")
tab_parent.add(tab3, text="    Create Fusion  ")


########WIDGETS for TAB0 #####################
Heading=ttk.Label(tab0,text="iGEM Tel Aviv 2019      TailOR-Swift",font=('Calibri', '12','bold'))
Heading.grid(row=0,column=0, padx=15, pady=10,columnspan=10)
#Heading=ttk.Label(tab0,text="TailOR Swift",font=('Calibri', '12','bold'))
#Heading.grid(row=0,column=2, padx=15, pady=10,columnspan=3)
im = Image.open('tailorswift.jpg')
ph = ImageTk.PhotoImage(im)
TailImage= ttk.Label(tab0,image=ph)
TailImage.grid(row=1,column=0,columnspan=10,rowspan=10)
InstructionLabel= ttkwidgets.LinkLabel(tab0, text="For software instruction press here",
link="https://2019.igem.org/Team:TAU_Israel",
normal_color=StColors.purple,
hover_color=StColors.torkiz,
clicked_color=StColors.light_grey)
InstructionLabel.grid(row=15,column=2, padx=15, pady=10,columnspan=4)




#tab_parent.add(tab0, text="profile", image=ph, compound=tk.TOP)

########WIDGETS for TAB1 #####################
Heading=ttk.Label(tab1,text="Tail Finder",font=('Calibri', '12','bold'))
Heading.grid(row=0,column=0, padx=15, pady=10,columnspan=2)

SearchBYLable=ttk.Label(tab1, text='Search By: ')
SearchBYEntry =ttk.Combobox(tab1, font=('Calibri', '12'), width=25)
#SearchBYEntry.bind('<<ComboboxSelected>>', lambda event:SearchBYEntry._click_combo())
SearchBYEntry['values'] = ('Host', 'Virus')
SearchBYLable.grid(row=1,column=0, padx=15, pady=10,columnspan=2, sticky=tk.W)
SearchBYEntry.grid(row=2,column=0, padx=15, pady=0,columnspan=2, sticky=tk.W+tk.E)

EmptyRow=ttk.Label(tab1)
EmptyRow.grid(row=3,column=0, padx=15, pady=0)


NameLabel=ttk.Label(tab1, text='Insert Name: ')
NameEntry=ttk.Entry(tab1, width=27,font=('Calibri', '12'))
NameLabel.grid(row=4,column=0, padx=15, pady=10,columnspan=2, sticky=tk.W, ipadx=100)
NameEntry.grid(row=5,column=0, padx=15, pady=0,columnspan=2, sticky=tk.W+tk.E)

EmptyRow=ttk.Label(tab1)
EmptyRow.grid(row=6,column=0, padx=15, pady=0)

explanationLabel=ttk.Label(tab1,text='Choose the required action')
explanationLabel.grid(row=7,column=0, columnspan=2,sticky=tk.W, padx=15, pady=10)
PerfectMatchButton=ttk.Button(tab1,text='Perfect Switch',width=12,command=perfect_match_finder_command)
PerfectMatchButton.grid(row=8,column=0)
FutionMatchButton=ttk.Button(tab1,text='Tails search',width=12,command=tail_searcher)
FutionMatchButton.grid(row=8,column=1)
DBsearchButton=ttk.Button(tab1,text='DB search',width=12,command=db_searcher)
DBsearchButton.grid(row=8,column=2)
style.configure("TButton", focuscolor=StColors.purple)

########WIDGETS for TAB2 #####################
Heading=ttk.Label(tab2,text="Find AA Sequence for Gene",font=('Calibri', '12','bold'))
Heading.grid(row=0,column=0, padx=15, pady=10,columnspan=10)

GeneCodeLable=ttk.Label(tab2, text='Gene Code: ')
GeneCodeEntry=ttk.Entry(tab2, width=27,font=('Calibri', '12'))
GeneCodeLable.grid(row=1,column=0, padx=15, pady=10,columnspan=2, sticky=tk.W)
GeneCodeEntry.grid(row=2,column=0, padx=15, pady=0,columnspan=2, sticky=tk.W+tk.E)

EmptyRow=ttk.Label(tab2)
EmptyRow.grid(row=3,column=0,columnspan=2, padx=15, pady=0)

AASeqButton=ttk.Button(tab2,text='Find',width=12,command=find_seq)
AASeqButton.grid(row=4,column=1)

########WIDGETS for TAB3 #####################


Heading=ttk.Label(tab3,text="Find the best location for fusion",font=('Calibri', '12','bold'))
Heading.grid(row=0,column=0, padx=15, pady=10,columnspan=10)

EmptyRow=ttk.Label(tab3)
EmptyRow.grid(row=1,column=0,columnspan=2, padx=15, pady=0)

FusionButton=ttk.Button(tab3,text='Create Fusion',width=12,command=create_fusion)
FusionButton.grid(row=2,column=5,columnspan=10)

tab_parent.pack(expand=1, fill='both')

gui.mainloop()


