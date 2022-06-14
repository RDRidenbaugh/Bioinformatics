import sys
import re
import random

#arguments inputed after calling the script in terminal
input = sys.argv[1] #your fasta formatted input file
output = sys.argv[2] #what you would like the the output file to be called
seed = sys.argv[3] #set seed for reproducibility. Must be an integer

random.seed(seed) #setting seed
infile = open(input, 'r')

#Dictonary of lists containing the iupac codings R, Y, S, W, K, M, B, D, H, and V
iupac = {'R': ['A', 'G'], 'Y': ['C', 'T'], 'S': ['C', 'G'], 'W': ['A', 'T'], 'K': ['G', 'T'], \
'M': ['A', 'C'], 'B': ['C', 'G', 'T'], 'D': ['A', 'G', 'T'], 'H': ['A', 'C', 'T'], 'V': ['A', 'C', 'G']}

with open(output, 'w') as f:
    for line in infile:
        if line.startswith('>'):
            f.write(line)
        else:
            #lines 23-48 replaces iupac codings with a random basepair based on the coding
            seq_list = list(line.upper())
            new_seq = ''
            for base in seq_list:
                if base == 'R':
                    new_seq = new_seq + random.choice(iupac['R'])
                elif base == 'Y':
                    new_seq = new_seq + random.choice(iupac['Y'])
                elif base == 'S':
                    new_seq = new_seq + random.choice(iupac['S'])
                elif base == 'W':
                    new_seq = new_seq + random.choice(iupac['W'])
                elif base == 'K':
                    new_seq = new_seq + random.choice(iupac['K'])
                elif base == 'M':
                    new_seq = new_seq + random.choice(iupac['M'])
                elif base == 'B':
                    new_seq = new_seq + random.choice(iupac['B'])
                elif base == 'D':
                    new_seq = new_seq + random.choice(iupac['D'])
                elif base == 'H':
                    new_seq = new_seq + random.choice(iupac['H'])
                elif base == 'V':
                    new_seq = new_seq + random.choice(iupac['V'])
                else:
                    new_seq = new_seq + base
            f.write(new_seq)
f.close()
infile.close()
