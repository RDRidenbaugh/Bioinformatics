#  Author: Ryan Ridenbaugh
#  python v.3.8.10
#  v.2.0

# REQUIRED FILE STRUCTURE
# . 					    Base Directory
# ├── /model                Directory containing blast search subject sequences in fasta format
# │     └──  query.fasta
# │ 
# ├── /non_model
# │     └──  subject.fasta        
# │ 
# ├── /blast_out       	    Directory containing *_out.txt blast output in outfmt6 format
# │     ├──  initial/       model -> non_model
# │     └──  reciprocal/    non_model -> model
# │ 
# └── reciprocal_blast.py 	

import os
import sys

argument = sys.argv[1]
base = sys.argv[2]

non_model_db = "makeblastdb -in "+base+"non_model/GCF_021901455.1_iyNeoLeco1.1_protein.fasta -dbtype prot -parse_seqids -out "+base+"non_model/db/nmodel_db" #EDIT CMD FOR YOUR PROJECT
model_db = "makeblastdb -in "+base+"model/FlyBase.fasta -dbtype prot -parse_seqids -out "+base+"model/db/model_db" #EDIT CMD FOR YOUR PROJECT
blastp = "blastp -query "+base+"model/FlyBase.fasta -db "+base+"non_model/db/nmodel_db -out "+base+"blast_out/initial/initial.txt -evalue 1e-10 -outfmt 6" #EDIT CMD FOR YOUR PROJECT

if argument == "-d":
    print("Generating databases..")
    os.system(non_model_db)
    os.system(model_db)
elif argument == "-dp":
    print("Generating databases and running BLAST+ for the initial search..")
    os.system(non_model_db)
    os.system(model_db)
    print("Running the initial BLAST+ search..")
    os.system(blastp)
elif argument == "-p":
    print("Running the initial BLAST+ search..")
    os.system(blastp)

parse = {}
for filename in os.listdir(base+"blast_out/initial"):
    nmodel_dict = {}
    with open(base+"non_model/GCF_021901455.1_iyNeoLeco1.1_protein.fasta", "r") as f: #EDIT PATH FOR YOUR PROJECT
        temp_list = []
        for line in f.readlines():
            if line.startswith('>'):
                temp_str = line.split()
                temp_list.append([temp_str[0], []])
            else: 
                temp_list[-1][1].append(line)
    nmodel_dict = {id: ''.join(seq).upper() for id, seq in temp_list}
    with open(os.path.join(base+"blast_out/initial/", filename), "r") as b_out_a:
        temp_list = []
        for line in b_out_a.readlines():
            temp_line = line.split()
            temp_list.append(temp_line)
for list in temp_list:
    if list[0] not in parse.keys():
        parse[list[0]] = [list[1]]
    elif list[0] in parse.keys():
        if list[1] not in parse[list[0]]:
            parse[list[0]].append(list[1])
for key in parse.keys():
    with open(base+"non_model/parsed/"+key+".fasta", "w") as parse_a:
        for id in parse[key]:
            parse_a.write('%s%s' % (">"+id+"\n", nmodel_dict[">"+id]))
for filename in os.listdir(base+"non_model/parsed/"):
        os.system("blastp -query "+base+"non_model/parsed/"+filename+" -db "+base+"model/db/model_db -out "+base+"blast_out/reciprocal/"+filename+".txt -evalue 1e-10 -outfmt 6") #EDIT CMD FOR YOUR PROJECT
