import os.path
import sys
import numpy as np

print("Starting to parse...")

def left(s, amount):
    return s[:amount]

def right(s, amount):
    return s[-amount:]

def mid(s, offset, amount):
    return s[offset:offset+amount]

fname = "dcon.out"
file_str = "./"+fname
if os.path.isfile(file_str):
    f = open(file_str, "r")
    lines = f.readlines()
    f.close()

str_search="Total Energy Eigenvalues:"
for i, line in enumerate(lines):
    line_clean = line.strip(" \t\n\r")
    if left(line_clean,len(str_search)) == str_search:
        i=i+4
        lsplit = lines[i].split()
        while len(lsplit) > 1:
            print(lsplit[1])
            i=i+1
            lsplit = lines[i].split()
