#!/usr/bin/python3
from rdkit import Chem
import numpy as np
from rdkit import DataStructs
from rdkit.Chem import rdMolDescriptors
import multiprocessing as mp

LIMIT=False
METALS=("Ti","Pd","Pt","Fe","Co","Ni","Cu","Zn","Ge","Po","Cd")
def matchBest(mainmols,subdb,top,useHyds=False,update=2.5e4):
    global LIMIT
    T=int(top)
    L=len(mainmols)
    scores=[[0 for _ in range(T)] for _ in mainmols]
    matches=[["" for _ in range(T)] for _ in mainmols]
    M=[max(sc) for sc in scores]
    m=[min(sc) for sc in scores]
    mainvecs=[Chem.RDKFingerprint(mainmol,useHs=useHyds) for mainmol in mainmols]
    checks=[False for _ in range(L)]
    K=0
    for s in subdb:
        compmol=Chem.MolFromSmiles(s)
        if not compmol:
            print("Omitting molecule:",s)
            continue
        compvec=Chem.RDKFingerprint(compmol,useHs=useHyds)
        closenesses=[DataStructs.TanimotoSimilarity(mainvec,compvec) for mainvec in mainvecs]
        #print(s,closeness,m)
        if all(checks): break
        for i in range(L):
            if checks[i]: continue
            if closenesses[i] > m[i]:
                for j in range(1,T):
                    if scores[i][j]<closenesses[i]:
                        scores[i][j-1]=scores[i][j]
                        matches[i][j-1]=matches[i][j]
                    else:
                        break
                else:
                    j+=1
                j-=1
            #print("Update:",j,"for",scores,closeness,"to get: ",end="")
            scores[i][j]=closenesses[i]
            matches[i][j]=s
            #print(scores)
            m[i]=scores[i][0]
            M[i]=scores[i][-1]
            if m[i]==1:  checks[i]=True
        K+=1
        if K%update==0:
            print(K,"molecules parsed.\t("+str(round((100*K)/len(subdb),2))+"% done)")
            reses=[list(zip(matches[i],scores[i])) for i in range(L)]
            for x in reses: print(x)
            print("------------------\n\n")
    return scores,matches
class MolDB:
    def __init__(self,infile,colID=1,refr=100000,dropmetals=False,**kwargs):
        self.smiles=[]
        curf=open(infile,"r")
        delim=' ' if "delim" not in kwargs else kwargs["delim"]
        K=0
        metskip=0
        for l in curf:
            l=l[0:len(l)-1] #Remove the '\n' in the end
            l=str.split(l,delim)
            if dropmetals:
                if any([metal in l[colID-1] for metal in METALS]):
                    metskip+=1
                    K+=1
                    continue
            self.smiles.append(l[colID-1])
            K+=1
            if K%refr==0: print(K,"molecules loaded")
        print(len(self.smiles),"molecules loaded into DB")
        if dropmetals: print(metskip,"SMILES skipped due to presence of heavy metals")


    def similarity(self,smis,top,useHyds=False,update=1e5):
        mainmols=[Chem.MolFromSmiles(smi) for smi in smis]
        return matchBest(mainmols,self.smiles,top,useHyds,update)
    def similarityMP(self,smis,top,procs=None,useHyds=False,update=1e5):
        global POOL,LIMIT
        mainmols=[Chem.MolFromSmiles(smi) for smi in smis]
        if not procs: procs=mp.cpu_count()
        frag=int(len(self.smiles)//procs)
        rem=len(self.smiles)%procs
        results=POOL.starmap(matchBest,[(mainmols,self.smiles[k:k+frag],top,useHyds) for k in range(0,len(self.smiles),frag)])
        if rem:
            print("Reminder correction:",rem)
            results.append(matchBest(mainmols,self.smiles[len(self.smiles)-rem:],top,useHyds))
        odata=results[0]
        for i in range(1,len(results)):
            for j in range(len(odata[0])):
                odata[0][j]+=results[i][0][j]
                odata[1][j]+=results[i][1][j]

        packed=list(zip(odata[0],odata[1]))
        packed=[tuple(zip(t[0],t[1])) for t in packed]
        #print("\nPacked:",packed)
        packed=[sorted(l,key=lambda x: x[0],reverse=True)[0:3] for l in packed]
        return packed
    '''
        T=int(top)
        scores=[0 for _ in range(T)]
        matches=["" for _ in range(T)]
        M=max(scores)
        m=min(scores)
        mainvec=Chem.RDKFingerprint(mainmol,useHs=useHyds)
        K=0
        for s in self.smiles:
            compmol=Chem.MolFromSmiles(s)
            if not compmol:
                print("Omitting molecule:",s)
                continue
            compvec=Chem.RDKFingerprint(compmol,useHs=useHyds)
            closeness=DataStructs.TanimotoSimilarity(mainvec,compvec)
            #print(s,closeness,m)
            if closeness > m:
                for j in range(1,T):
                    if scores[j]<closeness:
                        scores[j-1]=scores[j]
                        matches[j-1]=matches[j]
                    else:
                        break
                else:
                    j+=1
                j-=1
                #print("Update:",j,"for",scores,closeness,"to get: ",end="")
                scores[j]=closeness
                matches[j]=s
                #print(scores)
                m=scores[0]
                M=scores[-1]
                if m==1: break
            K+=1
            if K%update==0:
                print(K,"molecules parsed.\t("+str(round((100*K)/len(self.smiles),2))+"% done)")
                print(list(zip(matches,scores)))
        return scores,matches
    '''
#if __name__ == '__main__':
POOL = mp.Pool(mp.cpu_count())
