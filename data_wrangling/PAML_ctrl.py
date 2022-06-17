# python PAML maker

import os
from re import findall

A = open("PAML_A.job", "w")
B = open("PAML_B.job", "w")
C = open("PAML_C.job", "w")
D = open("PAML_D.job", "w")
E = open("PAML_E.job", "w")
F = open("PAML_F.job", "w")
G = open("PAML_G.job", "w")
H = open("PAML_H.job", "w")
I = open("PAML_I.job", "w")
J = open("PAML_J.job", "w")
free = open("PAML_free.job", "w")
prime = open("PAML_prime.sh", "w")

prime.write("""sbatch PAML_A.job
sbatch PAML_B.job
sbatch PAML_C.job
sbatch PAML_D.job
sbatch PAML_E.job
sbatch PAML_F.job
sbatch PAML_G.job
sbatch PAML_H.job
sbatch PAML_I.job
sbatch PAML_J.job
sbatch PAML_free.job""")
A.write("""#!/bin/sh
#SBATCH -t 14-00:00:00
#SBATCH --job-name=PAML_A
#SBATCH -n 1
#SBATCH -p normal
#SBATCH --account=coa_cli242_uksr
#SBATCH --mail-type ALL
#SBATCH --mail-user rdri228@g.uky.edu
#SBATCH --error=SLURM_JOB_%j.err
#SBATCH --output=SLURM_JOB_%j.out

module load Miniconda3
source /mnt/gpfs3_amd/share/apps/amd/Miniconda3/etc/profile.d/conda.sh
conda activate paml
 
""")
B.write("""#!/bin/sh
#SBATCH -t 14-00:00:00
#SBATCH --job-name=PAML_B
#SBATCH -n 1
#SBATCH -p normal
#SBATCH --account=coa_cli242_uksr
#SBATCH --mail-type ALL
#SBATCH --mail-user rdri228@g.uky.edu
#SBATCH --error=SLURM_JOB_%j.err
#SBATCH --output=SLURM_JOB_%j.out

module load Miniconda3
source /mnt/gpfs3_amd/share/apps/amd/Miniconda3/etc/profile.d/conda.sh
conda activate paml
 
""")
C.write("""#!/bin/sh
#SBATCH -t 14-00:00:00
#SBATCH --job-name=PAML_C
#SBATCH -n 1
#SBATCH -p normal
#SBATCH --account=coa_cli242_uksr
#SBATCH --mail-type ALL
#SBATCH --mail-user rdri228@g.uky.edu
#SBATCH --error=SLURM_JOB_%j.err
#SBATCH --output=SLURM_JOB_%j.out

module load Miniconda3
source /mnt/gpfs3_amd/share/apps/amd/Miniconda3/etc/profile.d/conda.sh
conda activate paml

""")
D.write("""#!/bin/sh
#SBATCH -t 14-00:00:00
#SBATCH --job-name=PAML_D
#SBATCH -n 1
#SBATCH -p normal
#SBATCH --account=coa_cli242_uksr
#SBATCH --mail-type ALL
#SBATCH --mail-user rdri228@g.uky.edu
#SBATCH --error=SLURM_JOB_%j.err
#SBATCH --output=SLURM_JOB_%j.out

module load Miniconda3
source /mnt/gpfs3_amd/share/apps/amd/Miniconda3/etc/profile.d/conda.sh
conda activate paml

""")
E.write("""#!/bin/sh
#SBATCH -t 14-00:00:00
#SBATCH --job-name=PAML_E
#SBATCH -n 1
#SBATCH -p normal
#SBATCH --account=coa_cli242_uksr
#SBATCH --mail-type ALL
#SBATCH --mail-user rdri228@g.uky.edu
#SBATCH --error=SLURM_JOB_%j.err
#SBATCH --output=SLURM_JOB_%j.out

module load Miniconda3
source /mnt/gpfs3_amd/share/apps/amd/Miniconda3/etc/profile.d/conda.sh
conda activate paml

""")
F.write("""#!/bin/sh
#SBATCH -t 14-00:00:00
#SBATCH --job-name=PAML_F
#SBATCH -n 1
#SBATCH -p normal
#SBATCH --account=coa_cli242_uksr
#SBATCH --mail-type ALL
#SBATCH --mail-user rdri228@g.uky.edu
#SBATCH --error=SLURM_JOB_%j.err
#SBATCH --output=SLURM_JOB_%j.out

module load Miniconda3
source /mnt/gpfs3_amd/share/apps/amd/Miniconda3/etc/profile.d/conda.sh
conda activate paml

""")
G.write("""#!/bin/sh
#SBATCH -t 14-00:00:00
#SBATCH --job-name=PAML_G
#SBATCH -n 1
#SBATCH -p normal
#SBATCH --account=coa_cli242_uksr
#SBATCH --mail-type ALL
#SBATCH --mail-user rdri228@g.uky.edu
#SBATCH --error=SLURM_JOB_%j.err
#SBATCH --output=SLURM_JOB_%j.out

module load Miniconda3
source /mnt/gpfs3_amd/share/apps/amd/Miniconda3/etc/profile.d/conda.sh
conda activate paml

""")
H.write("""#!/bin/sh
#SBATCH -t 14-00:00:00
#SBATCH --job-name=PAML_H
#SBATCH -n 1
#SBATCH -p normal
#SBATCH --account=coa_cli242_uksr
#SBATCH --mail-type ALL
#SBATCH --mail-user rdri228@g.uky.edu
#SBATCH --error=SLURM_JOB_%j.err
#SBATCH --output=SLURM_JOB_%j.out

module load Miniconda3
source /mnt/gpfs3_amd/share/apps/amd/Miniconda3/etc/profile.d/conda.sh
conda activate paml

""")
I.write("""#!/bin/sh
#SBATCH -t 14-00:00:00
#SBATCH --job-name=PAML_I
#SBATCH -n 1
#SBATCH -p normal
#SBATCH --account=coa_cli242_uksr
#SBATCH --mail-type ALL
#SBATCH --mail-user rdri228@g.uky.edu
#SBATCH --error=SLURM_JOB_%j.err
#SBATCH --output=SLURM_JOB_%j.out

module load Miniconda3
source /mnt/gpfs3_amd/share/apps/amd/Miniconda3/etc/profile.d/conda.sh
conda activate paml

""")
J.write("""#!/bin/sh
#SBATCH -t 14-00:00:00
#SBATCH --job-name=PAML_J
#SBATCH -n 1
#SBATCH -p normal
#SBATCH --account=coa_cli242_uksr
#SBATCH --mail-type ALL
#SBATCH --mail-user rdri228@g.uky.edu
#SBATCH --error=SLURM_JOB_%j.err
#SBATCH --output=SLURM_JOB_%j.out

module load Miniconda3
source /mnt/gpfs3_amd/share/apps/amd/Miniconda3/etc/profile.d/conda.sh
conda activate paml

""")
free.write("""#!/bin/sh
#SBATCH -t 14-00:00:00
#SBATCH --job-name=PAML_free
#SBATCH -n 1
#SBATCH -p normal
#SBATCH --account=coa_cli242_uksr
#SBATCH --mail-type ALL
#SBATCH --mail-user rdri228@g.uky.edu
#SBATCH --error=SLURM_JOB_%j.err
#SBATCH --output=SLURM_JOB_%j.out

module load Miniconda3
source /mnt/gpfs3_amd/share/apps/amd/Miniconda3/etc/profile.d/conda.sh
conda activate paml

""")
def PAML_ctl(transcripts, flavor, trees, control, out):
    ctl = """seqfile = {seq_file}
treefile = {tree_file}
outfile = {outfile}

noisy = 0
verbose = 0
runmode = 0

seqtype = 1
CodonFreq = 2
clock = 0
model = {model}

NSsites = 0 
icode = 0 

fix_kappa = 0 
kappa = 2 
fix_omega = {fix_omega} 
omega = 1

fix_alpha = 1
alpha = .0
Malpha = 0
ncatG = 4

getSE = 0
RateAncestor = 0
method = 0"""
    for file in os.listdir(transcripts):
        if findall(".fasta", file):
            if flavor == "A": # ωo = ωp = ωk
                XM = file.split(".")
                with open(control+XM[0]+"_A"+".ctl", "w") as A_ctl:
                    A_ctl.write(ctl.format(seq_file = transcripts+file, tree_file = trees+"5kb_unrooted_I.tre", outfile = out+XM[0]+"_A"+".txt", model = "0", fix_omega = "0"))
                A.write("codeml "+control+XM[0]+"_A"+".ctl\n")
            elif flavor == "B": # ωo = ωp, ωk
                XM = file.split(".")
                with open(control+XM[0]+"_B"+".ctl", "w") as B_ctl:
                    B_ctl.write(ctl.format(seq_file = transcripts+file, tree_file = trees+"5kb_unrooted_II.tre", outfile = out+XM[0]+"_B"+".txt", model = "2", fix_omega = "0"))
                B.write("codeml "+control+XM[0]+"_B"+".ctl\n")
            elif flavor == "C": # ωo = ωk, ωp
                XM = file.split(".")
                with open(control+XM[0]+"_C"+".ctl", "w") as C_ctl:
                    C_ctl.write(ctl.format(seq_file = transcripts+file, tree_file = trees+"5kb_unrooted_III.tre", outfile = out+XM[0]+"_C"+".txt", model = "2", fix_omega = "0"))
                C.write("codeml "+control+XM[0]+"_C"+".ctl\n")
            elif flavor == "D": # ωo, ωp = ωk
                XM = file.split(".")
                with open(control+XM[0]+"_D"+".ctl", "w") as D_ctl:
                    D_ctl.write(ctl.format(seq_file = transcripts+file, tree_file = trees+"5kb_unrooted_IV.tre", outfile = out+XM[0]+"_D"+".txt", model = "2", fix_omega = "0"))
                D.write("codeml "+control+XM[0]+"_D"+".ctl\n")
            elif flavor == "E": # ωo, ωp, ωk
                XM = file.split(".")
                with open(control+XM[0]+"_E"+".ctl", "w") as E_ctl:
                    E_ctl.write(ctl.format(seq_file = transcripts+file, tree_file = trees+"5kb_unrooted_V.tre", outfile = out+XM[0]+"_E"+".txt", model = "2", fix_omega = "0"))
                E.write("codeml "+control+XM[0]+"_E"+".ctl\n")
            elif flavor == "F":
                XM = file.split(".")
                with open(control+XM[0]+"_F"+".ctl", "w") as F_ctl:
                    F_ctl.write(ctl.format(seq_file = transcripts+file, tree_file = trees+"5kb_unrooted_II.tre", outfile = out+XM[0]+"_F"+".txt", model = "2", fix_omega = "1"))
                F.write("codeml "+control+XM[0]+"_F"+".ctl\n")
            elif flavor == "G":
                XM = file.split(".")
                with open(control+XM[0]+"_G"+".ctl", "w") as G_ctl:
                    G_ctl.write(ctl.format(seq_file = transcripts+file, tree_file = trees+"5kb_unrooted_III.tre", outfile = out+XM[0]+"_G"+".txt", model = "2", fix_omega = "1"))
                G.write("codeml "+control+XM[0]+"_G"+".ctl\n")
            elif flavor == "H":
                XM = file.split(".")
                with open(control+XM[0]+"_H"+".ctl", "w") as H_ctl:
                    H_ctl.write(ctl.format(seq_file = transcripts+file, tree_file = trees+"5kb_unrooted_IV.tre", outfile = out+XM[0]+"_H"+".txt", model = "2", fix_omega = "1"))
                H.write("codeml "+control+XM[0]+"_H"+".ctl\n")
            elif flavor == "I":
                XM = file.split(".")
                with open(control+XM[0]+"_I"+".ctl", "w") as I_ctl:
                    I_ctl.write(ctl.format(seq_file = transcripts+file, tree_file = trees+"5kb_unrooted_VI.tre", outfile = out+XM[0]+"_I"+".txt", model = "2", fix_omega = "1"))
                I.write("codeml "+control+XM[0]+"_I"+".ctl\n")
            elif flavor == "J":
                XM = file.split(".")
                with open(control+XM[0]+"_J"+".ctl", "w") as J_ctl:
                    J_ctl.write(ctl.format(seq_file = transcripts+file, tree_file = trees+"5kb_unrooted_VII.tre", outfile = out+XM[0]+"_J"+".txt", model = "2", fix_omega = "1"))
                J.write("codeml "+control+XM[0]+"_J"+".ctl\n")
            elif flavor == "free":
                XM = file.split(".")
                with open(control+XM[0]+"_free"+".ctl", "w") as free_ctl:
                    free_ctl.write(ctl.format(seq_file = transcripts+file, tree_file = trees+"5kb_unrooted_I.tre", outfile = out+XM[0]+"_free"+".txt", model = "1", fix_omega = "0"))
                free.write("codeml "+control+XM[0]+"_free"+".ctl\n")

PAML_ctl("../transcripts/", "A", "../trees/","../ctl/", "../paml_out/")
PAML_ctl("../transcripts/", "B", "../trees/","../ctl/", "../paml_out/")
PAML_ctl("../transcripts/", "C", "../trees/","../ctl/", "../paml_out/")
PAML_ctl("../transcripts/", "D", "../trees/","../ctl/", "../paml_out/")
PAML_ctl("../transcripts/", "E", "../trees/","../ctl/", "../paml_out/")
PAML_ctl("../transcripts/", "F", "../trees/","../ctl/", "../paml_out/")
PAML_ctl("../transcripts/", "G", "../trees/","../ctl/", "../paml_out/")
PAML_ctl("../transcripts/", "H", "../trees/","../ctl/", "../paml_out/")
PAML_ctl("../transcripts/", "I", "../trees/","../ctl/", "../paml_out/")
PAML_ctl("../transcripts/", "J", "../trees/","../ctl/", "../paml_out/")
PAML_ctl("../transcripts/", "free", "../trees/","../ctl/", "../paml_out/")

A.close()
B.close()
C.close()
D.close()
E.close()
F.close()
G.close()
H.close()
I.close()
J.close()
free.close()
prime.close()