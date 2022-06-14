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

# https://stackoverflow.com/questions/46921346/find-highest-depth-when-given-a-file-path
def get_depth(path, depth=0):
    if not os.path.isdir(path): return depth
    maxdepth = depth
    for entry in os.listdir(path):
        fullpath = os.path.join(path, entry)
        maxdepth = max(maxdepth, get_depth(fullpath, depth + 1))
    return maxdepth

import sys
from re import findall
from re import sub

iter_count = 0
depth = get_depth(sys.argv[1])
while iter_count <= depth+1:
    for dirpath, dirnames, filenames in os.walk(sys.argv[1], topdown= True):
        iter_count = iter_count + 1
        if dirnames or filenames:
            for str in dirnames:
                if findall(r"[\/*\-\ \,\?\/!]+", str):
                    new_path = sub(r"[\/*\-\ \,\?\/!]+", "_", str)
                    os.rename(os.path.join(dirpath, str), os.path.join(dirpath, new_path))
            for str in filenames:
                if findall(r"[\/*\-\ \,\?\/!]+", str):
                    new_path = sub(r"[\/*\-\ \,\?\/!]+", "_", str)
                    os.rename(os.path.join(dirpath, str), os.path.join(dirpath, new_path))
