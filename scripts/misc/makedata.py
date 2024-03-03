#!/usr/bin/env python3
# Generate data for molecules by counting elements and their occurences

import sys
from rdkit import Chem

fname=sys.argv[1]
f=open(fname,"r")
elems=('C','c','N','n','O','P','S','F','Cl','Br','I')
err=0
for l in f:
    l=l[0:len(l)-1]
    sec1=str.split(l,'\t')
    molname=sec1[1]
    smiles=sec1[0]
    print(l)
    m=Chem.MolFromSmiles(smiles)
    if not m:
        err+=1
        continue
    fp=Chem.RDKFingerprint(m)
    fp=[str(i) for i in fp]
    cs=[]
    for e in elems:
        cs.append(str(smiles.count(e)))
    print(molname,smiles,''.join(fp),' '.join(cs))
print("Completed with "+str(err)+" errors")
