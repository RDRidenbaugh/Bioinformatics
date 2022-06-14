# v.0

import os

gff = []
XMid = {}
genome = {}
XMloc_neg = {}
XMloc_pos = {}
scaffold_key = {}
converter = {"NC_060260.1":"Chromosome1","NC_060261.1":"Chromosome2","NC_060262.1":"Chromosome3","NC_060263.1":"Chromosome4","NC_060264.1":"Chromosome5","NC_060265.1":"Chromosome6","NC_060266.1":"Chromosome7"}
def transcipt_loc(XM, gff):
    with open (XM, "r") as f:
        for line in f.readlines():
            temp_split = line.replace("\n", "").split("\t")
            XMid[temp_split[0]] = temp_split[1]
    f.close()
    with open(gff, "r") as f:
        for line in f.readlines():
            if line.startswith("NC"):
                temp_split  = line.split()
                if temp_split[2] == "exon":
                    meta = temp_split[8].split(";")
                    meta = meta[1].split("-")
                    for key in XMid.keys():
                        if meta[-1] == key:
                            scaffold_key[meta[-1]] = converter[temp_split[0]]
                            if temp_split[6] == "-":
                                if meta[-1] not in XMloc_neg.keys():
                                    XMloc_neg[meta[-1]] = [[temp_split[3], temp_split[4]]]
                                else:
                                    XMloc_neg[meta[-1]].insert(0, [temp_split[3], temp_split[4]])
                            else:
                                if meta[-1] not in XMloc_pos.keys():
                                    XMloc_pos[meta[-1]] = [[temp_split[3], temp_split[4]]]
                                else:
                                    XMloc_pos[meta[-1]].append([temp_split[3], temp_split[4]])
    f.close()
def transcript_extract(path,out): 
    for file in os.listdir(path):
        temp_species = file.split('_')[0]
        with open(path+file, 'r') as f:
            temp_list = []
            for line in f.readlines():
                if line.startswith('>'):
                    temp_index = line.find(" ")
                    temp_fasta = line[0:temp_index].replace(">","")
                    temp_list.append([temp_fasta, []])
                else: 
                    temp_list[-1][1].append(line)
            temp_dict = {id: ''.join(seq).upper().replace('\n', '') for id, seq in temp_list}
        genome[temp_species] = temp_dict
        f.close()
    for k in XMid.keys():
        f = open(out+k+".fasta", "w")
        if k in XMloc_pos.keys():
            for species in genome.keys():
                temp_str = ""
                for list in XMloc_pos[k]:
                    start = int(list[0])
                    stop = int(list[-1])
                    temp_substring = genome[species][scaffold_key[k]][start-1:stop]
                    temp_str = temp_str + temp_substring.strip()
                f.write(">"+k+"_"+species+" "+XMid[k]+"\n"+temp_str+"\n")
        elif k in XMloc_neg.keys():
            for species in genome.keys():
                temp_str = ""
                for list in XMloc_neg[k]:
                    start = int(list[0])
                    stop = int(list[-1])
                    temp_substring = genome[species][scaffold_key[k]][start-1:stop]
                    temp_str = temp_str + temp_substring.strip()
                    reverse = ""
                for base in reversed(temp_str):
                    if base == "A":
                        reverse = reverse + "T"
                    elif base == "T":
                        reverse = reverse + "A"
                    elif base == "G":
                        reverse = reverse + "C"
                    elif base == "C":
                        reverse = reverse + "G"
                    else:
                        reverse = reverse + base
                f.write(">"+k+"_"+species+" "+XMid[k]+"\n"+reverse+"\n")
        f.close()
        
transcipt_loc("../transcripts/mechano_transcripts.txt", "../GCF_021901455.1_iyNeoLeco1.1_genomic.gff")
transcript_extract("../reference/", "../transcripts/")
# print(XMid.keys())
# print(genome.keys())
