#!/usr/bin/python3
# Count number of each atom type class. Each class can contain multilpe atom types and the counts are grouped by class

import os
import sys

if len(sys.argv):
    files=[]
    fp=open(sys.argv[1],"r")
    for l in fp:
        print(l)
        l=str.split(l.replace("\t"," "),' ')[0]
        files.append(l)
else:
    files=os.listdir()
print(files)
registers=[["CT1","CT2","CT3","CT3x"],["CT"],["C3"],["CT2x","CT1x"],["CP2","CTL3"],["HA","HB","HP"],["HE1","HE2"],["HAL3","HAL2","HAL1","HL"],["H","HC"],["OH1","OH"],["OHL"],["NH1","NH2","N"],["NH3","NP"],["CPT","CA","CY"],["CN"],["NC"],["CAS"],["NX","NY"],["CAG"],["NC2","NC2T"],["CP1","CP3"],["CF1","CF2"],["CF3"],["FA"],["F1","F2","F3"],["CLAL"],["CLAR"],["BRAR"],["BRAL"],["I"],["CXT1","CXT2"],["CTL2","CTL1"],["PL"],["SL"],["CE1A","CE1B"],["CE2A","CE2B"],["C","CC","CD","CL"],["O","OC"],["OCL"],["OB"],["OD"],["OS"],["OSL"],["O2L"],["HF1","HF2"],["HOL"],["SM","S"],["HS","HD"],["SX"]]
finvecs=[]
pdbfiles=[]
for f in files:
    if not f.endswith(".pdb"): continue
    pdbfiles.append(f)
    print(f)
    curf=open(f,"r")
    tempv=[0 for _ in range(len(registers))]
    for l in curf:
        l=l[0:len(l)-1] #Remove the '\0' in the end
        if not (l.startswith("ATOM") or l.startswith("HETATM")): continue
        atype=str.strip(l[12:17])
        ind=-1
        for i in range(len(registers)):
            if atype in registers[i]:
                ind=i
                break
        if ind==-1:
            print("Atom type: '"+atype+"' not found in register")
            continue
        tempv[ind]+=1
    finvecs.append(tempv)
for i in range(len(finvecs)):
    print(pdbfiles[i],*finvecs[i])
if len(finvecs)!=len(pdbfiles):
    print("Warn: Vector count doesn't match file count")
    print(len(finvecs),len(pdbfiles))
