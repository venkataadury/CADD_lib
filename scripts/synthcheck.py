#!/usr/bin/python3
from rdkit import Chem,DataStructs
from syba.syba import SybaClassifier
import matplotlib.pyplot as plt
import sys
import numpy as np
from math import sqrt,cos,sin,acos,asin,exp,log,pi

sqr=lambda x: x*x
distance=lambda p1,p2: sqrt(sqr(p2[0]-p1[0])+sqr(p2[1]-p1[1]))

def loadMolDataFromSMILESFile(fileobj):
    mols=[]
    fnames=[]
    print("Loading mols from SMILES file... ")
    for f in fileobj:
        f=f[:len(f)-1]; #Remove the '\n' in the end
        f=f.split()
        print(f)
        fnames.append(f[1])
        mols.append(Chem.MolFromSmiles(f[0],sanitize=True))
    print("Done")

    fails=sum([1 for i in mols if i is None])
    print(len(mols),"molecules loaded with",fails,"failures")
    print("Cleaning ... ",end="")
    mols=[mol for mol in mols if mol]
    print("done")
    return mols,fnames

def makeSimilarityMatrix(set1,set2,symmetric=False,ignoreonecopy=False): #Set1 and Set2 are Fingerprint Lists
    N1=len(set1); N2=len(set2)
    simmat=np.zeros((N1,N2))
    if symmetric and N1!=N2: raise ValueError("Symmetric matrix requested, but set sizes are not same!")
    for i in range(N1):
        copy=0
        if symmetric: k=i+1
        else: k=0
        for j in range(k,N2):
            simmat[i][j]=DataStructs.TanimotoSimilarity(set1[i],set2[j])
            if symmetric: simmat[j][i]=simmat[i][j]
            else:
                if simmat[i][j]>0.99 and copy<1:
                    copy+=1
                    simmat[i][j]=0
        if (i+1)%100==0: print(i,"molecules processed of",N1)
    return simmat

fn=sys.argv[1]
fl=open(fn,"r");
mols,fnames=loadMolDataFromSMILESFile(fl)
syba = SybaClassifier()
syba.fitDefaultScore()
succ=0
for i,name in enumerate(fnames):
    acc=syba.predict(Chem.MolToSmiles(mols[i]))>0
    print(name,acc)
    if acc: succ+=1

print("Synthesizable Fraction:",(succ/len(mols)))
