# v.1.0

import os
from re import findall
from scipy import stats

def PAML_parse(path, XMkey, out):
    dnr = set()
    results = {}
    for file in os.listdir(path):
        XM = ""
        if os.stat(path+file).st_size == 0:
            temp_split = file.split(".")
            temp_split = temp_split[0].split("_")
            dnr.add(temp_split[0]+"_"+temp_split[1])
        else:
            temp_split = file.split(".")
            temp_split = temp_split[0].split("_")
            XM = temp_split[0]+"_"+temp_split[1]
            flavour = temp_split[-1]
            with open(path+file, "r") as f:
                for line in f.readlines():
                    if findall(r"ln[lL]", line):
                        temp_split = line.split()
                        temp_lnL = temp_split[4]
                        temp_np = temp_split[3][0:-2] # Works if np = two didgets
                        if XM not in results.keys():
                            results[XM] = {flavour : [temp_lnL, temp_np]}
                        elif flavour not in results[XM].keys():
                             results[XM][flavour] = [temp_lnL, temp_np]
                        else:
                            results[XM][flavour].append([temp_lnL, temp_np])
            f.close()
    with open(out, "w") as o:
        o.write("XM\tgene_id\tA-free\tA-B\tA-C\tA-D\tB-E\tC-E\tB-F\tC-G\tD-H\tE-I\tE-J\n")
        with open (XMkey, "r") as f:
            XMid = {}
            for line in f.readlines():
                temp_split = line.replace("\n", "").split("\t") 
                XMid[temp_split[0]] = temp_split[1]
        f.close()
        for key, value in XMid.items():
            if key[0:-2] in dnr:
                pass
            else:
                Afree_p = stats.chi2.sf(2*(float(results[key[0:-2]]["free"][0]) - float(results[key[0:-2]]["A"][0])), int(results[key[0:-2]]["free"][-1]) - int(results[key[0:-2]]["A"][-1]))
                AB_p = stats.chi2.sf(2*(float(results[key[0:-2]]["B"][0]) - float(results[key[0:-2]]["A"][0])), int(results[key[0:-2]]["B"][-1]) - int(results[key[0:-2]]["A"][-1]))
                AC_p = stats.chi2.sf(2*(float(results[key[0:-2]]["C"][0]) - float(results[key[0:-2]]["A"][0])), int(results[key[0:-2]]["C"][-1]) - int(results[key[0:-2]]["A"][-1]))
                AD_p = stats.chi2.sf(2*(float(results[key[0:-2]]["D"][0]) - float(results[key[0:-2]]["A"][0])), int(results[key[0:-2]]["D"][-1]) - int(results[key[0:-2]]["A"][-1]))
                BE_p = stats.chi2.sf(2*(float(results[key[0:-2]]["E"][0]) - float(results[key[0:-2]]["B"][0])), int(results[key[0:-2]]["E"][-1]) - int(results[key[0:-2]]["B"][-1]))
                CE_p = stats.chi2.sf(2*(float(results[key[0:-2]]["E"][0]) - float(results[key[0:-2]]["C"][0])), int(results[key[0:-2]]["E"][-1]) - int(results[key[0:-2]]["C"][-1]))
                BF_p = stats.chi2.sf(2*(float(results[key[0:-2]]["B"][0]) - float(results[key[0:-2]]["F"][0])), int(results[key[0:-2]]["B"][-1]) - int(results[key[0:-2]]["F"][-1]))
                CG_p = stats.chi2.sf(2*(float(results[key[0:-2]]["C"][0]) - float(results[key[0:-2]]["G"][0])), int(results[key[0:-2]]["C"][-1]) - int(results[key[0:-2]]["G"][-1]))
                DH_p = stats.chi2.sf(2*(float(results[key[0:-2]]["D"][0]) - float(results[key[0:-2]]["H"][0])), int(results[key[0:-2]]["D"][-1]) - int(results[key[0:-2]]["H"][-1]))
                EI_p = stats.chi2.sf(2*(float(results[key[0:-2]]["E"][0]) - float(results[key[0:-2]]["I"][0])), int(results[key[0:-2]]["E"][-1]) - int(results[key[0:-2]]["I"][-1]))
                EJ_p = stats.chi2.sf(2*(float(results[key[0:-2]]["E"][0]) - float(results[key[0:-2]]["J"][0])), int(results[key[0:-2]]["E"][-1]) - int(results[key[0:-2]]["J"][-1]))
                o.write(key+"\t"+value+"\t"+str(Afree_p)+"\t"+str(AB_p)+"\t"+str(AC_p)+"\t"+str(AD_p)+"\t"+str(BE_p)+"\t"+str(CE_p)+"\t"+str(BF_p)+"\t"+str(CG_p)+"\t"+str(DH_p)+"\t"+str(EI_p)+"\t"+str(EJ_p)+"\n")
    o.close()
                        
PAML_parse("../paml_out/", "../transcripts/mechano_transcripts.txt", "../lnL_results.tab")