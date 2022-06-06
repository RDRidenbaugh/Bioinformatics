# v.3.2

gff = []
def gff_read(path):
    with open(path, "r") as f:
        for line in f.readlines():
            if line.startswith("NC"):
                temp_split  = line.split()
                if temp_split[2] == "CDS":
                    gff.append(line.split("\t"))

qtl = []
genome = {}
peak_genes = {}
metadata = {}
c_key = {"1":"NC_060260.1", "2":"NC_060261.1", "3":"NC_060262.1", "4":"NC_060263.1", "5":"NC_060264.1", "6":"NC_060265.1", "7":"NC_060266.1"}
def compare_range(start1, end1, start2, end2):
    return max(max((end2-start1), 0) - max((end2-end1), 0) - max((start2-start1), 0), 0)
def qtl_extract(path_a, path_b, path_c):
    with open(path_a, "r") as f:
        for line in f.readlines():
            temp_split = line.split()
            start = temp_split[5].split("_")
            stop = temp_split[6].split("_")
            qtl.append([temp_split[0]+"_"+temp_split[1],temp_split[3], start[1], stop[1]])
    with open(path_b, 'r') as f:
        temp_list = []
        for line in f.readlines():
            if line.startswith('>'):
                temp_index = line.find(" ")
                temp_fasta = line[0:temp_index].replace(">","")
                temp_list.append([temp_fasta, []])
            else: 
                temp_list[-1][1].append(line)
        genome = {id: ''.join(seq).upper().replace('\n', '') for id, seq in temp_list}
    for list_a in qtl:
        for list_b in gff:
            d = list_b[8].split(";")
            if c_key[list_a[1]] == list_b[0]:
                overlap = compare_range(int(list_a[2]), int(list_a[3]), int(list_b[3]), int(list_b[4]))
                if overlap != 0:
                    if list_a[0] not in peak_genes.keys():
                        peak_genes[list_a[0]] = set()
                        temp_split = list_b[8].split(";")
                        temp_split = temp_split[0].split("-")
                        peak_genes[list_a[0]] = set()
                        peak_genes[list_a[0]].add(temp_split[-1])
                        metadata[temp_split[-1]] = ["[Last_Exon_Range:"+list_b[3]+"-"+list_b[4]+" "+d[1]+" "+d[5]+" "+d[6]+"]"+"\n"]
                    else:
                        temp_split = list_b[8].split(";")
                        temp_split = temp_split[0].split("-")
                        peak_genes[list_a[0]].add(temp_split[-1])
                        metadata[temp_split[-1]] = ["[Last_Exon_Range:"+list_b[3]+"-"+list_b[4]+" "+d[1]+" "+d[5]+" "+d[6]+"]"+"\n"]
    output = open(path_c, "w")
    for key in peak_genes.keys():
        for transcript in peak_genes[key]:
            output.write(">"+transcript+" "+key+" "+metadata[transcript][0]+genome[transcript]+"\n")
    output.close()

gff_read("GCF_021901455.1_iyNeoLeco1.1_genomic.gff")
qtl_extract("color_qtl_loc.txt", "GCF_021901455.1_iyNeoLeco1.1_protein.faa", "color_qtl_genes.faa")