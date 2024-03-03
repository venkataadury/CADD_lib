#!/usr/bin/env python3
import numpy as np
from rdkit import Chem
import sys
from syba import syba

infile=open(sys.argv[1],"r")
passfile=open(sys.argv[2] if len(sys.argv)>2 else "success.smi","w")
errfile=open(sys.argv[3] if len(sys.argv)>3 else "error.smi","w")
failfile=open(sys.argv[4] if len(sys.argv)>4 else "failures.smi","w")

sc=syba.SybaClassifier()
sc.fitDefaultScore()
print("Fitted Score",flush=True)
def validationFunction(mol):
    return sc.predict(mol)>0

for l in infile:
    l=l[0:len(l)-1].split() # Remove '\n' in the end
    sm=l[0]
    name = l[1] if len(l)>1 else "Molecule_unnamed.pdb"
    try:
        res=validationFunction(sm)
        if res: passfile.write(sm+" "+name+"\n")
        else: failfile.write(sm+" "+name+"\n")
    except:
        errfile.write(sm+" "+name+"\n")
passfile.close()
errfile.close()
failfile.close()

