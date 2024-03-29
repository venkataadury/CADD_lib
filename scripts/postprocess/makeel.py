#!/usr/bin/python3
#
# Replace all atom names with element names in a PDB file

import sys
import os

file_names = [fn for fn in os.listdir(".") if fn.endswith(".pdb")] #Select all PDB files in this folder
for fl in file_names:
    newfnam=str.split(fl,'.')[0]+"_el.pdb"
    outf=open(newfnam,"w")
    inf=open(fl,"r")
    print(fl)
    for ln in inf:
        ln=ln[0:len(ln)-1] #Remove the '\n' at the end of variable:-ln
        if ln[0:4]=="ATOM" or ln[0:6]=="HETATM": outf.write(ln[0:12]+ln[len(ln)-4:]+ln[16:]+"\n")
        else: outf.write(ln+"\n")
    outf.close()
    inf.close()
