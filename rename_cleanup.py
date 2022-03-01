# v.1.0
# python v.3.8.10

# AUTHOR: 
# Ryan D. Ridenbaugh

# DESCRIPTION:
# Python script that cleans up directory and file names by recursively moving though \
# directories from the top down replacing spaces with "_". The script can target additional symbols \
# by modifying the regular expression in the findall() function.

# REQUIRED ARGUMENTS:
# Path to the top directory is required after calling the script - "python3 rename_cleanup.py /PATH"

import os
import sys
from re import findall

for i, (dirpath, dirnames, filenames) in enumerate(os.walk(sys.argv[1])):
    if i == 0:   
        if dirnames:
            for str in dirnames:
                if findall(r" ", str):
                    old_path = dirpath+str
                    new_path = dirpath+str.replace(" ", "_")
                    os.renames(old_path, new_path)
        if filenames:
            for str in filenames:
                if findall(r" ", str):
                    old_path = dirpath+str
                    print(old_path)
                    new_path = dirpath+str.replace(" ", "_")
                    os.renames(old_path, new_path)
    elif i > 0:
        if dirnames:
            for str in dirnames:
                if findall(r" ", str):
                    old_path = dirpath+"/"+str
                    new_path = dirpath+"/"+str.replace(" ", "_")
                    os.renames(old_path, new_path)
        if filenames:
            for str in filenames:
                if findall(r" ", str):
                    old_path = dirpath+"/"+str
                    print(old_path)
                    new_path = dirpath+"/"+str.replace(" ", "_")
                    os.renames(old_path, new_path)
