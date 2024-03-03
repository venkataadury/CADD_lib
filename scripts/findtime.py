#!/usr/bin/python3

import sys
from math import sqrt

fl=open(sys.argv[1],"r")
p1=float(sys.argv[2])
p2=float(sys.argv[3])
mD1,mD2=1e10,1e10
mD=[1e10]*5
t1,t2=0,0
t=[0]*5
for l in fl:
    l=str.split(l.strip(),' ')
    d1,d2=abs(float(l[1])-p1),abs(float(l[2])-p2)
    d3=sqrt(d1**2+d2**2)
    if d1<mD1:
        mD1=d1
        t1=float(l[0])
        continue
    if d2<mD2:
        mD2=d2
        t2=float(l[0])
        continue
    if d3<mD[-1]:
        ind=0
        for i in range(len(mD)-1,-1,-1):
            if mD[i]<=d3:
                ind=i
                break
        mD=mD[0:ind]+[d3]+mD[ind:len(mD)-1]
        t=t[0:ind]+[float(l[0])]+t[ind:len(mD)-1]
print(t1,t2)
print(mD)
print(t)
