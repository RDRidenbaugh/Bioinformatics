gff = []
def gff_read(path):
    with open(path, "r") as f:
        for line in f.readlines():
            if line.startswith("NC"):
                temp_split  = line.split()
                if temp_split[2] == "gene":
                    gff.append(line.split())

qtl = []
scaffolds = {}
c_key = {"1":"NC_060260.1", "2":"NC_060261.1", "3":"NC_060262.1", "4":"NC_060263.1", "5":"NC_060264.1", "6":"NC_060265.1", "7":"NC_060266.1"}
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
        scaffolds = {id: ''.join(seq).upper().replace('\n', '') for id, seq in temp_list}
    output = open(path_c, "w")
    for list_a in qtl:
        for list_b in gff:
            if c_key[list_a[1]] == list_b[0]:
                range_a = range(int(list_a[2]), int(list_a[3]))
                range_b = range(int(list_b[3]), int(list_b[4]))
                overlap = [i for i in range_a if i in range_b]
                if overlap:
                    start = int(list_b[3])-1
                    stop = int(list_b[4])-1
                    temp_substring = scaffolds[list_b[0]][start:stop]
                    output.write(">"+list_a[0]+" "+start+"-"+stop+" "+list_b[-1]+"\n"+temp_substring+"\n")
                else:
                    pass    
    output.close()

gff_read("../non_model/GCF_021901455.1_iyNeoLeco1.1_genomic.gff")
qtl_extract("../non_model/color_qtl_loc.txt", "../non_model/GCF_021901455.1_iyNeoLeco1.1_genomic.fna", "../non_model/color_qtl_genes.fasta")