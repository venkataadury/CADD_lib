#!/usr/bin/python3

import sys
valences=[[],["HA","HB","HE1","HE2","HC","HP","HD","H","HOL","HS","HL","HF1","HF2","OC","OB","O","OD","O2L","OCL","F1","F2","F3","FA","HAL3","HAL2","HAL1","CLAL","CLAR","BRAL","BRAR","I","NC"],["OH1","OHL","NX","CN","OS","OSL","NC2T","S","SX","SM"],["CA","CPT","CY","CAS","CD","CE1A","CE2A","CE1B","CE2B","C","CC","NH1","NH2","NH3","N","NP","NY","NC2","CAG","CL"],["CT1","CT2","CT3","CT3x","CF1","CF2","CF3","CT","CP2","CP1","CP3","CXT1","CXT2","C3","SL","PL","CT1x","CT2x","CTL3","CTL2","CTL1"]]
allv=[]; 
for el in valences: allv+=el
def readPDB(f):
    atomnums=dict()
    bonds=dict()
    of=open(f,"r")
    for l in of:
        if (l.startswith("ATOM") or l.startswith("HETATM")): 
            atomnums[int(l[8:11])]=str.strip(l[12:16])
        elif l.startswith("CONECT"):
            n=((len(str.strip(l))-6)//5)-1
            atn=int(l[7:11])
            bonds[atn]=n
            if atomnums[atn] in valences[n]: continue
            else:
                if atomnums[atn] in allv: print("WARN: Atom number '"+str(atn)+"' of atom-type",atomnums[atn],"has unsatsified valence:",n)
                else: print("NOTE: Atom-type '"+atomnums[atn]+"' is not found in database.")
        else: continue
    return atomnums,bonds
atomnums,bonds=readPDB(sys.argv[1])
print("Done with",sys.argv[1])
#for el in atomnums.keys(): print(el,atomnums[el]+": ",bonds[el])
