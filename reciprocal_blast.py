# REQUIRED FILE STRUCTURE
# . 					
# ├── /model                Directory containing blast search subject sequences in fasta format
# │     └──  query.fasta
# ├── /non_model
# │     └──  subject.fasta /db        
# ├── /blast_out       	    Directory containing *_out.txt blast output in outfmt6 format
# └── reciprocal_blast.jl 	
#
# OPTIONAL FILES

import os
import sys
import subprocess
from re import findall

argument = sys.argv[1]
base = os.path.dirname(sys.argv[2])+"/"
if argument == "-i": #interactive local mode
    blast_db_a = input("Do you need to build a BLAST database for the non-model organism? (y/n): ")
# while True:
    if blast_db_a.lower() == "y" or blast_db_a.lower() == "yes":
        print("yes")
        in_file_a = input("Name of the file for database creation: ")
        in_type_a = input("What type of file is this? (nucl/prot): ")
        nmodel_db = input("What would you like to name your blast database: ")
        try:
            subprocess.run(["makeblastdb", "-in", base+"non_model/"+in_file_a, "-dbtype", in_type_a, "-parse_seqids", "-out", base+"non_model/db/"+nmodel_db])
        except FileNotFoundError as error:
            print("An error has occured. Please check that NCBI BLAST+ is installed correctly and the information entered is accurate.")
            print(error)
    else:
        print("no")
        pass
    blast_db_b = input("Do you need to build a BLAST database for the model organism? (y/n): ")
    if blast_db_b.lower() == "y" or blast_db_b.lower() == "yes":
        in_file_b = input("Name of the file for database creation: ")
        in_type_b = input("What type of file is this? (nucl/prot): ")
        model_db = input("What would you like to name your blast database: ")
        try:
            subprocess.run(["makeblastdb", "-in", base+"model/"+in_file_b, "-dbtype", in_type_b, "-parse_seqids", "-out", base+"model/db/"+model_db])
        except FileNotFoundError as error:
            print("An error has occured. Please check that NCBI BLAST+ is installed correctly and the information entered is accurate.")
            print(error)
    flavor_a = input("What type of BLAST are you performing?: ")
    model_query_a = input("Name of your starting query file: ")
    if 'nmodel_db' not in locals():
        nmodel_db = input("Name of your starting BLAST database: ")
    evalue_a = input("evalue cutoff: ")
    out_a = input("What would you like to name your outfile?: ")
    try:
        subprocess.run([flavor_a,"-query", base+"model/"+model_query_a, "-db", base+"non_model/db/"+nmodel_db, "-out", base+"blast_out/"+out_a+".txt", "-evalue", evalue_a, "-outfmt", "6"])
    except FileNotFoundError as error:
        print("An error has occured. Please check that NCBI BLAST+ is installed correctly and the information entered is accurate.")
        print(error)
else: #remote  cluster mode
    control = []
    with open(base+"control.txt", "r") as ctrl:
        for line in ctrl.readlines():
            control.append(line.split())
    ctrl.close()
    if argument == "-d":
        print("Generating databases..")
        subprocess.run(control[0])
        subprocess.run(control[1])
    elif argument == "-dp":
        print("Generating databases and running BLAST+ for the initial search..")
        subprocess.run(control[0])
        subprocess.run(control[1])
        subprocess.run(control[2])
    elif argument == "-p":
        print("Running the initial BLAST+ search..")
        subprocess.run(control[2])
if 'out_a' not in locals():
    parse = {}
    for filename in os.listdir(base+"blast_out/initial"):
        nmodel_dict = {}
        with open(base+"non_model/test_nmodel.fasta", "r") as f: 
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

else:       
    print("Parsing functionality for interactive local mode coming soon")
print("Running the reciprocal BLAST+ search..")
for filename in os.listdir(base+"non_model/parsed/"):
        subprocess.run(["blastp", "-query", base+"non_model/parsed/"+filename, "-db", base+"model/db/test_model_db", "-out", base+"blast_out/reciprocal/"+filename+".txt", "-evalue", "1e-10", "-outfmt", "6"])

