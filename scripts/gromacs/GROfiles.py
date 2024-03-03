#!/usr/bin/python3
import numpy as np
from math import sin,cos,asin,acos,pi,sqrt

def PBCWrap(dV,pbcx,pbcy=None,pbcz=None):
    if pbcy is None: pbcy=pbcx
    if pbcz is None: pbcz=pbcy
    if pbcx:
        dV[0]%=pbcx
        if 2*dV[0]>=pbcx: dV[0]-=pbcx
    if pbcy:
        dV[1]%=pbcy
        if 2*dV[1]>=pbcy: dV[1]-=pbcy
    if pbcz:
        dV[2]%=pbcz
        if 2*dV[2]>=pbcz: dV[2]-=pbcz
    return dV
def LJPotential(sigma,epsilon,r):
    s=sigma/r
    return 4*epsilon*(s**12-s**6)
def CoulombPotential(q1,q2,r):
    return (138.935*q1*q2)/r
def AtomicPotential(a1,a2):
    d=a1.getDistanceTo(a2)
    try:
        return LJPotential((a1.sigma+a2.sigma)/2,sqrt(a1.epsilon*a2.epsilon),d)+CoulombPotential(a1.charge,a2.charge,d)
    except:
        print("Atom names:",a1.name,a2.name)
        raise ValueError("SIGMA and EPSILON not available for these atoms")

class Atom:
    def __init__(self,name,x=0,y=0,z=0,resnum=1,resname="UNK"):
        self.name=name.strip()
        self.pos=np.array([x,y,z])
        self.resnum=resnum
        self.resname=resname.strip()

    def getX(self): return self.pos[0]
    def getY(self): return self.pos[1]
    def getZ(self): return self.pos[2]

    def setX(self,x): self.pos[0]=x
    def setY(self,y): self.pos[1]=y
    def setZ(self,z): self.pos[2]=z

    def getCharge(self): return self.charge
    def setCharge(self,q): self.charge=q
    def getSigma(self): return self.sigma
    def setSigma(self,s): self.sigma=s
    def getEpsilon(self): return self.epsilon
    def setEpsilon(self,e): self.epsilon=e

    def getPositionVector(self): return self.pos

    def setResidueNumber(self,n): self.resnum=n
    def getResidueNumber(self): return self.resnum
    def setResidue(self,res): self.resname=res
    def getResidue(self): return self.resname

    def getVectorTo(self,targ): return targ.pos-self.pos
    def getVectorFrom(self,src): return self.pos-src.pos
    def getDistanceTo(self,sec): return np.linalg.norm(self.pos-sec.pos)
    def getDistanceFrom(self,sec): return np.linalg.norm(self.pos-sec.pos)

    def getVectorToPBC(self,targ,pbcx,pbcy=None,pbcz=None):
        dV=self.getVectorTo(targ)
        return PBCWrap(dV,pbcx,pbcy,pbcz)
    def getVectorFromPBC(self,targ,pbcx,pbcy=None,pbcz=None):
        dV=self.getVectorFrom(targ)
        return PBCWrap(dV,pbcx,pbcy,pbcz)
    def getDistanceToPBC(self,sec,pbcx,pbcy=None,pbcz=None): return np.linalg.norm(self.getVectorFromPBC(sec,pbcx,pbcy,pbcz))
    def getDistanceFromPBC(self,sec,pbcx,pbcy=None,pbcz=None): return np.linalg.norm(self.getVectorFromPBC(sec,pbcx,pbcy,pbcz))
        

class GROFile:
    #upd is the frequency to write loaded frame count 
    #reduction is the set of functions that operate on the data to reduce it if the file is too large. They can act only on one frame at a time, and the entire frame atom-data is wiped out from memory after all reduction operations are completed unless any of the functions themselves preserves a portion of it.
    def __init__(self,fnam,upd=250,reduction=None):
        self.frames=[]
        self.frametitles=[]
        atoms=[]
        newfile=open(fnam,"r")
        title=True
        count=False
        dump=False
        N=0
        K=0
        for l in newfile:
            if dump:
                title=True
                dump=False
                continue
            if title:
                self.frametitles.append(l[0:len(l)-1])
                title=False
                count=True
                continue
            if count:
                N=int(l.strip())
                K=0
                count=False
                continue
            resnum=int(l[0:5])
            resname=l[5:10]
            atname=l[10:15]
            xc=float(l[22:28].strip())
            yc=float(l[30:36].strip())
            zc=float(l[37:44].strip())
            #print(xc,yc,zc)
            atoms.append(Atom(atname,xc,yc,zc,resnum,resname))
            K+=1
            if K>=N:
                dump=True
                if reduction is not None and reduction: atoms=[f(atoms) for f in reduction]
                self.frames.append(atoms)
                atoms=[]
                if upd and len(self.frames)%upd==0: print(len(self.frames),"frames loaded ...")
        print("Loaded file completely")
        print(len(self.frames),"frames with",len(self.frames[0]),"atoms per frame")

    def writeGROFile(self,filename):
        nf=open(filename,"w")
        for i in range(len(self.frames)):
            nf.write(self.frametitles[i]+"\n")
            nf.write(" "+str(len(self.frames[i]))+"\n")
            for j,atom in enumerate(self.frames[i]):
                nf.write("{0:5d}{1:<4s}{2:>5s} {3:5d}{4:>8s}{5: >8s}{6: >8s}".format(atom.resnum,atom.resname,atom.name,j,str(round(atom.getX(),3)),str(round(atom.getY(),3)),str(round(atom.getZ(),3)))+"\n")
            nf.write(" 0.000 0.000 0.000\n")
        nf.close()

    def addAtom(self,a,framenums=None):
        if framenums is None: framenums=list(range(len(self.frames)))
        for fi in framenums:
            self.frames[fi].append(a)

    #Apply the set of functions (frame-by-frame) to each frame and give the resultant list as output. Original data is preserved
    def produceCondensed(self,funcs):
        out=[]
        for f in funcs:
            fout=[]
            for frame in self.frames:
                fout.append(f(frame))
            out.append(fout)
        return out

import re
#Pre-written analysis code templates for full file data (non-reduced). Some of these templates can be used as reduction functions also
def condCompile(condS,condT):
    if condT is not None:
        try: condT=list(condT)
        except: condT=[condT]
    if condS is not None:
        try: condS=list(condS)
        except: condS=[condS]
    return condS,condT

#This function will find the closest target to each source atom (each atom matching source regexp) in ONE frame
PBCall=7.11 #Can be modified
def proximal(frame,src_type,targ_type,cut=0,condT=None,condS=None):
    global PBCall
    condS,condT=condCompile(condS,condT)
    dists=np.zeros((len(frame),len(frame)))
    results=[]
    acc,acc2=True,True
    for i in range(len(frame)):
        if condS:
            for cond in condS:
                if not cond(frame[i]):
                    acc=False
                    break
            else:
                acc=True
        if acc and re.match(src_type,frame[i].name):
            for j in range(len(frame)):
                if condT:
                    for cond in condT:
                        if not cond(frame[j],frame[i]):
                            acc2=False
                            break
                    else:
                        acc2=True
                if acc2 and re.match(targ_type,frame[j].name):
                    dists[i][j]=frame[i].getDistanceFromPBC(frame[j],PBCall)
                    if cut and dists[i][j]>cut: dists[i][j]=0
                    dists[j][i]=dists[i][j]
            mD=1e8
            mI=-1
            for k,d in enumerate(dists[i]):
                if k==i: continue
                if d and mD>d:
                    mD=d
                    mI=k
            if mI!=-1: results.append((i,mI,mD))
    return results

def angleBetween(a1,a2,a3):
    v1=a2.getVectorTo(a1)
    v2=a2.getVectorTo(a3)
    d=np.dot(v1,v2)
    #print("Info:",d,np.linalg.norm(v1),np.linalg.norm(v2))
    try: return acos(d/(np.linalg.norm(v1)*np.linalg.norm(v2)))
    except: return pi if d<0 else 0

def atom_position_extractor(frame):
    poses=[]
    for atom in frame:
        poses.append([atom.getPositionVector()])
    return poses


def cluster(frame,threshold=0.7,pointred=atom_position_extractor):
    global PBCall #Modify to set PBC boundaries
    points=pointred(frame)
    '''
    dists=[[threshold+1 for _ in range(len(points))] for _ in range(len(points))]
    for i in range(len(points)):
        for j in range(i,len(points)):
            if i==j:
                dists[i][j]=0
                continue
            else:
                dists[i][j]=np.linalg.norm(PBCWrap(points[i]-points[j],PBCall))
                dists[j][i]=dists[i][j]
    '''
    clusters=[[(i,points[i])] for i in range(len(points))]
    succ=True
    while succ:
        succ=False
        for i in range(len(clusters)):
            for j in range(i+1,len(clusters)):
                for n1 in clusters[i]:
                    for n2 in clusters[j]:
                        for p1 in n1[1]:
                            for p2 in n2[1]:
                                #if dists[n1[0]][n2[0]]<=threshold:
                                if np.linalg.norm(PBCWrap(p1-p2,PBCall))<=threshold:
                                    clusters=[clusters[k] for k in range(len(clusters)) if (k!=i and k!=j)]+[clusters[i]+clusters[j]]
                                    succ=True
                                    break
                            if succ: break
                        if succ: break
                    if succ: break
                if succ: break
            if succ: break
    return clusters


def assignParameters(frame,params):
    for atom in frame:
        if atom.name in params:
            atom.sigma=params[atom.name][0]
            atom.epsilon=params[atom.name][1]
            atom.charge=params[atom.name][2]
