from re import split, findall, sub
from time import strftime, localtime

print('Initializing function "name"', strftime("%Y-%m-%d %I:%M:%S %p", localtime()))
def name(protein_names):
    with open(protein_names, "r") as pn:
        temp_dict = {}
        for line in pn.readlines():
            line = line.strip("\n")
            if line.startswith("#"):
                continue
            else:
                temp_split = line.split("\t")
                total_split = line.replace("\t", " ").split(" ")
                t = len(total_split)-1
                temp_value = []
                # if t >= n:
                #     continue
                # else:
                if temp_split[1] == '':
                    temp_value = [["NA"]]
                elif len(temp_split[1]) == 14:
                    temp_value = [[temp_split[1]]]
                else:
                    temp_value = [temp_split[1].split(" ")]
                if temp_split[2] == '':
                    temp_value.append(["NA"])
                elif len(temp_split[2]) == 14:
                    temp_value.append([temp_split[2]])
                else:
                    temp_value.append(temp_split[2].split(" "))
                if temp_split[3] == '':
                    temp_value.append(["NA"])
                elif len(temp_split[3]) == 14:
                    temp_value.append([temp_split[3]])
                else:
                    temp_value.append(temp_split[3].split(" "))
                if temp_split[4] == '':
                    temp_value.append(["NA"])
                elif len(temp_split[4]) == 14:
                    temp_value.append([temp_split[4]])
                else:
                    temp_value.append(temp_split[4].split(" "))
                if temp_split[5] == '':
                    temp_value.append(["NA"])
                elif len(temp_split[5]) == 14:
                    temp_value.append([temp_split[5]])
                else:
                    temp_value.append(temp_split[5].split(" "))
                temp_dict[temp_split[0]] = temp_value
    pn.close()
    return temp_dict

print('Running function "name" for table_OGs_protein_names.txt', strftime("%Y-%m-%d %I:%M:%S %p", localtime()))
OG_names = name("../broccoli_gene_orthology/run1/dir_step3/table_OGs_protein_names.txt")
print('Function "name" executed!', strftime("%Y-%m-%d %I:%M:%S %p", localtime()))

print('Initializing function "build_db_ft"', strftime("%Y-%m-%d %I:%M:%S %p", localtime()))
def build_db_ft(path, feature):
    temp_dict = {}
    with open(path, "r") as f:
        for line in f.readlines():
            if line.startswith(feature):
                if feature == "mRNA":
                    temp_split = line.strip().split("\t")
                    temp_dict[temp_split[10]] = (temp_split[14], temp_split[10], temp_split[12], int(temp_split[18]), int(temp_split[17]), temp_split[5], temp_split[13])
                if feature == "CDS":
                    temp_split = line.strip().split("\t")
                    temp_dict[temp_split[10]] = (temp_split[14], temp_split[12], temp_split[10], int(temp_split[18]), int(temp_split[17]), temp_split[5], temp_split[13])
    f.close()
    return temp_dict

print('Running function "build_db_ft" for the five genomes', strftime("%Y-%m-%d %I:%M:%S %p", localtime()))
pine_db_ft = build_db_ft("../feature_table/GCF_021155775.1_iyNeoPine1.1_feature_table.txt", "CDS")
leco_db_ft = build_db_ft("../feature_table/GCF_021901455.1_iyNeoLeco1.1_feature_table.txt", "CDS")
virg_db_ft = build_db_ft("../feature_table/GCF_021901495.1_iyNeoVirg1.1_feature_table.txt", "CDS")
fabr_db_ft = build_db_ft("../feature_table/GCF_021155785.1_iyNeoFabr1.1_feature_table.txt", "CDS")
simi_db_ft = build_db_ft("../feature_table/GCF_021155765.1_iyDipSimi1.1_feature_table.txt", "CDS")
print('Function "build_db_ft" executed!', strftime("%Y-%m-%d %I:%M:%S %p", localtime()))

print('Initializing function "build_db_exon"', strftime("%Y-%m-%d %I:%M:%S %p", localtime()))
def build_db_exon(path, feature):
    XM_dict = {}
    with open(path, "r") as f:
        for line in f.readlines():
            if findall(feature, line) and not findall("partial=true", line):
                temp_split = line.strip("\n").split("\t")
                meta_split = split(r'[;=]',temp_split[-1])
                if meta_split[-1] not in XM_dict:
                    XM_dict[meta_split[-1]] = [[temp_split[3], temp_split[4]]]
                else:
                    XM_dict[meta_split[-1]].append([temp_split[3], temp_split[4]])
    f.close()
    for key in XM_dict.keys():
        exon_list = XM_dict[key]
        exon_len = len(exon_list)
        if exon_len == 1:
            temp_list_exon = [abs(int(exon_list[0][1])-int(exon_list[0][0]))+1]
            XM_dict[key] = temp_list_exon
        else:
            temp_list_exon = []
            temp_list_intron = []
            for sublist in exon_list:
                temp_len = abs(int(sublist[1])-int(sublist[0]))+1
                temp_list_exon.append(temp_len)
            for i, sublist in enumerate(exon_list):
                if i == 0:
                    temp_start = sublist[1]
                elif i < exon_len-1 and i != 0:
                    temp_len = abs(int(temp_start)-int(sublist[0]))
                    temp_list_intron.append(temp_len)
                    temp_start = sublist[1]
                else:
                    temp_len = abs(int(temp_start)-int(sublist[0]))
                    temp_list_intron.append(temp_len)
            XM_dict[key] = [temp_list_exon, temp_list_intron]
    # print(XM_dict["XM_046611274.1"]) [[19, 454, 118, 196, 148, 240, 120, 273, 257, 1134], [9831, 3487, 1246, 455, 806, 1067, 477, 598, 1637]]
    # print(XM_dict["XM_046611268.1"]) [[379, 454, 118, 196, 148, 240, 120, 273, 257, 1134], [6734, 3487, 1246, 455, 806, 1067, 477, 598, 1637]]
    return XM_dict

print('Running function "build_db_exon" for the five genomes', strftime("%Y-%m-%d %I:%M:%S %p", localtime()))
pine_db_exon = build_db_exon("../genome_gff/GCF_021155775.1_iyNeoPine1.1_genomic.gff", "exon")
leco_db_exon = build_db_exon("../genome_gff/GCF_021901455.1_iyNeoLeco1.1_genomic.gff", "exon")
virg_db_exon = build_db_exon("../genome_gff/GCF_021901495.1_iyNeoVirg1.1_genomic.gff", "exon")
fabr_db_exon = build_db_exon("../genome_gff/GCF_021155785.1_iyNeoFabr1.1_genomic.gff", "exon")
simi_db_exon = build_db_exon("../genome_gff/GCF_021155765.1_iyDipSimi1.1_genomic.gff", "exon")
print('Function "build_db_exon" executed!', strftime("%Y-%m-%d %I:%M:%S %p", localtime()))

print('Initializing function "longest"', strftime("%Y-%m-%d %I:%M:%S %p", localtime()))
def longest(list, species):
    if species == "pine":
        high = 0
        XP = ""
        for i in range(len(list)):
            if pine_db_ft[list[i]][3] > high:
                high = pine_db_ft[list[i]][3]
                XP = list[i]
    elif species == "leco":
        high = 0
        XP = ""
        for i in range(len(list)):
            if leco_db_ft[list[i]][3] > high:
                high = leco_db_ft[list[i]][3]
                XP = list[i]
    elif species == "virg":
        high = 0
        XP = ""
        for i in range(len(list)):
            if virg_db_ft[list[i]][3] > high:
                high = virg_db_ft[list[i]][3]
                XP = list[i]
    elif species == "fabr":
        high = 0
        XP = ""
        for i in range(len(list)):
            if fabr_db_ft[list[i]][3] > high:
                high = fabr_db_ft[list[i]][3]
                XP = list[i]
    elif species == "simi":
        high = 0
        XP = ""
        for i in range(len(list)):
            if simi_db_ft[list[i]][3] > high:
                high = simi_db_ft[list[i]][3]
                XP = list[i]
    return XP

print('Initializing function "OG_filter"', strftime("%Y-%m-%d %I:%M:%S %p", localtime()))
def OG_filter(dict):
    temp_dict = {}
    for key, list in dict.items():
        if len(list[0]) == 1 and len(list[1]) == 1 and len(list[2]) == 1 and len(list[3]) == 1 and len(list[4]) == 1:
            temp_dict[key] = (list[0][0], list[1][0], list[2][0], list[3][0], list[4][0])
        else:
            if len(list[0]) == 1:
                if list[0][0] == "NA":
                    temp_dict[key] = ("NA",)
                else:
                    temp_dict[key] = (list[0][0],)
            elif len(list[0]) > 1:
                temp_dict[key] = (longest(list[0], "pine"),)  
            if len(list[1]) == 1:
                if list[1][0] == "NA":
                    temp_dict[key] = temp_dict[key]+("NA",)
                else:
                    temp_dict[key] = temp_dict[key]+(list[1][0],)
            elif len(list[1]) > 1:
                temp_dict[key] = temp_dict[key]+(longest(list[1], "leco"),)
            if len(list[2]) == 1:
                if list[2][0] == "NA":
                    temp_dict[key] = temp_dict[key]+("NA",)
                else:
                    temp_dict[key] = temp_dict[key]+(list[2][0],)
            elif len(list[2]) > 1:
                temp_dict[key] = temp_dict[key]+(longest(list[2], "virg"),)
            if len(list[3]) == 1:
                if list[3][0] == "NA":
                    temp_dict[key] = temp_dict[key]+("NA",)
                else:
                    temp_dict[key] = temp_dict[key]+(list[3][0],)
            elif len(list[3]) > 1:
                temp_dict[key] = temp_dict[key]+(longest(list[3], "fabr"),)
            if len(list[4]) == 1:
                if list[4][0] == "NA":
                    temp_dict[key] = temp_dict[key]+("NA",)
                else:
                    temp_dict[key] = temp_dict[key]+(list[4][0],)
            elif len(list[4]) > 1:
                temp_dict[key] = temp_dict[key]+(longest(list[4], "simi"),)
    return temp_dict

print('Running function "OG_filter" for OG_names', strftime("%Y-%m-%d %I:%M:%S %p", localtime()))
OG_names = OG_filter(OG_names)
print('Function "OG_filter" executed!', strftime("%Y-%m-%d %I:%M:%S %p", localtime()))

print('Initializing function "pairwise"', strftime("%Y-%m-%d %I:%M:%S %p", localtime()))
def pairwise(key, dif):
    for i, o in enumerate(key):
        if i == 0:
            pine = int(pine_db_ft[o][3])
        elif i == 1:
            leco = int(leco_db_ft[o][3])
        elif i == 2:
            virg = int(virg_db_ft[o][3])
        elif i == 3:
            fabr = int(fabr_db_ft[o][3])
        elif i == 4:
            simi = int(simi_db_ft[o][3])
    if abs(pine - leco) >= dif: #or abs(pine - virg) >= n or abs(pine - fabr) >= n or abs(leco - virg) >= n or abs(leco - fabr) >= n or abs(virg - fabr) >= n:
        return True
    else:
        return False

print('Initializing function "OG_lp"', strftime("%Y-%m-%d %I:%M:%S %p", localtime()))
def OG_lp(dict, dif):
    temp_dict_a = {}
    temp_dict_b = {}
    for key in dict.keys():
        if dict[key][0] == "NA" or dict[key][1] == "NA" or dict[key][2] == "NA" or dict[key][3] == "NA" or dict[key][4] == "NA":
            temp_dict_a[key] = dict[key]
        elif pairwise(dict[key], dif):
            temp_dict_b[key] = dict[key]
        else: 
            continue
    return temp_dict_a, temp_dict_b

print('Running function "OG_lp" for OG_names', strftime("%Y-%m-%d %I:%M:%S %p", localtime()))
output = OG_lp(OG_names, 100)
print('Function "OG_lp" executed!', strftime("%Y-%m-%d %I:%M:%S %p", localtime()))

print('Initializing function "lp_output"', strftime("%Y-%m-%d %I:%M:%S %p", localtime()))
def lp_output(dict, path):
    with open(path, "w") as f:
        f.write("ortho_group\tLOC\tCHR\tpinetum_xm\tlecontei_xm\tvirginianus_xm\tfabricii_xm\tsimilis_xm\tpinetum_len\tlecontei_len\tvirginianus_len\tfabricii_len\tsimilis_len\tproduct\n")
        for key in dict.keys():
                f.write(key+"\t"+pine_db_ft[dict[key][0]][0]+"\t"+pine_db_ft[dict[key][0]][5]+"\t"+pine_db_ft[dict[key][0]][1]+"\t"+leco_db_ft[dict[key][1]][1]+"\t"+virg_db_ft[dict[key][2]][1]+"\t"+ fabr_db_ft[dict[key][3]][1]+"\t"+simi_db_ft[dict[key][4]][1]+"\t")
                f.write(str(pine_db_ft[dict[key][0]][3])+"\t"+str(leco_db_ft[dict[key][1]][3])+"\t"+str(virg_db_ft[dict[key][2]][3])+"\t"+str(fabr_db_ft[dict[key][3]][3])+"\t"+str(simi_db_ft[dict[key][4]][3])+"\t"+pine_db_ft[dict[key][0]][6]+"\n")
    f.close()

print('Initializing function "pa_output"', strftime("%Y-%m-%d %I:%M:%S %p", localtime()))
def pa_output(dict, path):
    with open(path, "w") as f:
        f.write("ortho_group\tpinetum_xp\tlecontei_xp\tvirginianus_xp\tfabricii_xp\tsimilis_xp\tchromosome\tproduct\n")
        for key in dict.keys():
            f.write(key+"\t")
            if dict[key][0] != "NA":
                f.write(pine_db_ft[dict[key][0]][1]+"\t")
            else:
                f.write("NA\t")
            if dict[key][1] != "NA":
                f.write(leco_db_ft[dict[key][1]][1]+"\t")
            else:
                f.write("NA\t")
            if dict[key][2] != "NA":
                f.write(virg_db_ft[dict[key][2]][1]+"\t")
            else:
                f.write("NA\t")
            if dict[key][3] != "NA":
                f.write(fabr_db_ft[dict[key][3]][1]+"\t")
            else:
                f.write("NA\t")
            if dict[key][4] != "NA":
                f.write(simi_db_ft[dict[key][4]][1]+"\t")
            else:
                f.write("NA\t")
            if dict[key][0] != "NA":
                f.write(pine_db_ft[dict[key][0]][5]+"\t")
                f.write(pine_db_ft[dict[key][0]][6]+"\n")
            elif dict[key][1] != "NA":
                f.write(leco_db_ft[dict[key][1]][5]+"\t")
                f.write(leco_db_ft[dict[key][1]][6]+"\n")
            elif dict[key][2] != "NA":
                f.write(virg_db_ft[dict[key][2]][5]+"\t")
                f.write(virg_db_ft[dict[key][2]][6]+"\n")
            elif dict[key][3] != "NA":
                f.write(fabr_db_ft[dict[key][3]][5]+"\t")
                f.write(fabr_db_ft[dict[key][3]][6]+"\n")
            elif dict[key][4] != "NA":
                f.write(simi_db_ft[dict[key][4]][6]+"\t")
                f.write(simi_db_ft[dict[key][4]][6]+"\n")
        f.close()

print('Outputting data using function "lp_output" and "pa_output"', strftime("%Y-%m-%d %I:%M:%S %p", localtime()))
lp_output(output[1], "../length_polymorphism.txt")
pa_output(output[0], "../presence_polymorphism.txt")
print('Function "lp_output" and "pa_output" executed!', strftime("%Y-%m-%d %I:%M:%S %p", localtime()))