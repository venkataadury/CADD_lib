#!/usr/bin/python3
# Calculate running averages from a file (assumed colums are 't' and 'x')

#Running averages
import numpy as np
import sys
DIVI=1
if len(sys.argv)>1: fname=sys.argv[1]
else: fname="anal.spc"
rawdata=np.loadtxt(fname)
ts=[]
ds=[]
for i in range(len(rawdata)):
    ts.append(rawdata[i][0]) #Timestamp
    ds.append(sum(rawdata[i][1:])/DIVI)
NFRAME=3500
ra=[]
nns=[]
rs=sum([ds[k] for k in range(NFRAME-1)])
for i in range(len(ds)-NFRAME):
    rs+=ds[NFRAME+i-1]
    nns.append(ts[i])
    ra.append(rs/NFRAME)
    rs-=ds[i]
outf=open("running_avg.spc","w")
for i in range(len(ra)):
    outf.write(str(nns[i])+" "+str(ra[i])+"\n")
outf.close()
