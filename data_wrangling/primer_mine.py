from time import strftime, localtime

print('Initializing function "name"', strftime("%Y-%m-%d %I:%M:%S %p", localtime()))
def name(protein_names, n):
    with open(protein_names, "r") as pn:
        temp_dict = {}
        for line in pn.readlines():
            line = line.strip("\n")
            if line.startswith("#"):
                pass
            else:
                temp_split = line.split("\t")
                total_split = line.replace("\t", " ").split(" ")
                t = len(total_split)-1
                temp_value = []
                if t >= n:
                    continue
                else:
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
OG_names = name("../broccoli_gene_orthology/run1/dir_step3/table_OGs_protein_names.txt", 15)
print('Function "name" executed!', strftime("%Y-%m-%d %I:%M:%S %p", localtime()))

print('Initializing function "build_db"', strftime("%Y-%m-%d %I:%M:%S %p", localtime()))
def build_db(path):
    temp_dict = {}
    with open(path, "r") as f:
        for line in f.readlines():
            if line.startswith("CDS"):
                temp_split = line.strip().split("\t")
                temp_dict[temp_split[10]] = (temp_split[14], temp_split[12], temp_split[10], int(temp_split[18]), int(temp_split[17]), temp_split[5], temp_split[13])
    f.close()
    return temp_dict

print('Running function "build_db" for the five genomes', strftime("%Y-%m-%d %I:%M:%S %p", localtime()))
pine_db = build_db("../feature_table/GCF_021155775.1_iyNeoPine1.1_feature_table.txt")
leco_db = build_db("../feature_table/GCF_021901455.1_iyNeoLeco1.1_feature_table.txt")
virg_db = build_db("../feature_table/GCF_021901495.1_iyNeoVirg1.1_feature_table.txt")
fabr_db = build_db("../feature_table/GCF_021155785.1_iyNeoFabr1.1_feature_table.txt")
simi_db = build_db("../feature_table/GCF_021155765.1_iyDipSimi1.1_feature_table.txt")
print('Function "build_db" executed!', strftime("%Y-%m-%d %I:%M:%S %p", localtime()))

print('Initializing function "longest"', strftime("%Y-%m-%d %I:%M:%S %p", localtime()))
def longest(list, species):
    if species == "pine":
        high = 0
        XP = ""
        for i in range(len(list)):
            if pine_db[list[i]][3] > high:
                high = pine_db[list[i]][3]
                XP = list[i]
    elif species == "leco":
        high = 0
        XP = ""
        for i in range(len(list)):
            if leco_db[list[i]][3] > high:
                high = leco_db[list[i]][3]
                XP = list[i]
    elif species == "virg":
        high = 0
        XP = ""
        for i in range(len(list)):
            if virg_db[list[i]][3] > high:
                high = virg_db[list[i]][3]
                XP = list[i]
    elif species == "fabr":
        high = 0
        XP = ""
        for i in range(len(list)):
            if fabr_db[list[i]][3] > high:
                high = fabr_db[list[i]][3]
                XP = list[i]
    elif species == "simi":
        high = 0
        XP = ""
        for i in range(len(list)):
            if simi_db[list[i]][3] > high:
                high = simi_db[list[i]][3]
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
def pairwise(key, n):
    for i, o in enumerate(key):
        if i == 0:
            pine = int(pine_db[o][3])
        elif i == 1:
            leco = int(leco_db[o][3])
        elif i == 2:
            virg = int(virg_db[o][3])
        elif i == 3:
            fabr = int(fabr_db[o][3])
        elif i == 4:
            simi = int(simi_db[o][3])
    if abs(pine - leco) >= n or abs(pine - virg) >= n or abs(pine - fabr) >= n or abs(leco - virg) >= n or abs(leco - fabr) >= n or abs(virg - fabr) >= n:
        return True
    else:
        return False

print('Initializing function "OG_lp"', strftime("%Y-%m-%d %I:%M:%S %p", localtime()))
def OG_lp(dict, n):
    temp_dict_a = {}
    temp_dict_b = {}
    for key in dict.keys():
        if dict[key][0] == "NA" or dict[key][1] == "NA" or dict[key][2] == "NA" or dict[key][3] == "NA" or dict[key][4] == "NA":
            temp_dict_a[key] = dict[key]
        elif pairwise(dict[key], n):
            temp_dict_b[key] = dict[key]
        else: 
            continue
    return temp_dict_a, temp_dict_b

print('Running function "OG_lp" for OG_names', strftime("%Y-%m-%d %I:%M:%S %p", localtime()))
output = OG_lp(OG_names, 67)
print('Function "OG_lp" executed!', strftime("%Y-%m-%d %I:%M:%S %p", localtime()))

print('Initializing function "lp_output"', strftime("%Y-%m-%d %I:%M:%S %p", localtime()))
def lp_output(dict, path):
    with open(path, "w") as f:
        f.write("ortho_group\tpinetum_xm\tlecontei_xm\tvirginianus_xm\tfabricii_xm\tsimilis_xm\tpinetum_len\tlecontei_len\tvirginianus_len\tfabricii_len\tsimilis_len\tproduct\n")
        for key in dict.keys():
                f.write(key+"\t"+pine_db[dict[key][0]][1]+"\t"+leco_db[dict[key][1]][1]+"\t"+virg_db[dict[key][2]][1]+"\t"+ fabr_db[dict[key][3]][1]+"\t"+simi_db[dict[key][4]][1]+"\t")
                f.write(str(pine_db[dict[key][0]][3])+"\t"+str(leco_db[dict[key][1]][3])+"\t"+str(virg_db[dict[key][2]][3])+"\t"+str(fabr_db[dict[key][3]][3])+"\t"+str(simi_db[dict[key][4]][3])+"\t"+pine_db[dict[key][0]][6]+"\n")
    f.close()

print('Initializing function "pa_output"', strftime("%Y-%m-%d %I:%M:%S %p", localtime()))
def pa_output(dict, path):
    with open(path, "w") as f:
        f.write("ortho_group\tpinetum_xp\tlecontei_xp\tvirginianus_xp\tfabricii_xp\tsimilis_xp\tproduct\n")
        for key in dict.keys():
            f.write(key+"\t")
            if dict[key][0] != "NA":
                f.write(pine_db[dict[key][0]][1]+"\t")
            else:
                f.write("NA\t")
            if dict[key][1] != "NA":
                f.write(leco_db[dict[key][1]][1]+"\t")
            else:
                f.write("NA\t")
            if dict[key][2] != "NA":
                f.write(virg_db[dict[key][2]][1]+"\t")
            else:
                f.write("NA\t")
            if dict[key][3] != "NA":
                f.write(fabr_db[dict[key][3]][1]+"\t")
            else:
                f.write("NA\t")
            if dict[key][4] != "NA":
                f.write(simi_db[dict[key][4]][1]+"\t")
            else:
                f.write("NA\t")
            if dict[key][0] != "NA":
                f.write(pine_db[dict[key][0]][6]+"\n")
            elif dict[key][1] != "NA":
                f.write(leco_db[dict[key][1]][6]+"\n")
            elif dict[key][2] != "NA":
                f.write(virg_db[dict[key][2]][6]+"\n")
            elif dict[key][3] != "NA":
                f.write(fabr_db[dict[key][3]][6]+"\n")
            elif dict[key][4] != "NA":
                f.write(simi_db[dict[key][4]][6]+"\n")
        f.close()

print('Outputting data using function "lp_output" and "pa_output"', strftime("%Y-%m-%d %I:%M:%S %p", localtime()))
lp_output(output[1], "../length_polymorphism.txt")
pa_output(output[0], "../presence_polymorphism.txt")
print('Function "lp_output" and "pa_output" executed!', strftime("%Y-%m-%d %I:%M:%S %p", localtime()))