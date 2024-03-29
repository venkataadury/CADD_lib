#!/usr/bin/python3
#
# (Mathematically) edit a given column in a file

import os
import sys
import argparse
from math import sin,cos,tan,asin,acos,atan,sqrt

parser=argparse.ArgumentParser(prog="coledit.py")
parser.add_argument(dest="infile",action="store",help="The input file")
parser.add_argument('-d',dest='delim',metavar='--delimiter',action="store",help="Column Delimiter",default=" ")
parser.add_argument('-c',dest='col',metavar='--column',action="store",help="Column Number. Can specify character range also (eg. -c25-27)",default=None)
parser.add_argument('-m',dest='mult',metavar='--multiplier',action="store",help="Multiplier",default=1)
parser.add_argument('-M',dest='math',metavar='--math',action="store",help="Any mathematical expression (If you do not want to use plain multiplier). Use 'x' as the data value. (eg. 3*x+5)",default=None)
parser.add_argument('--multiplyfirst', dest='mfirst', action='store_const',const="1", default="0",help='If using multiplier and expression together, perform multiplication first, then apply the expression')
parser.add_argument('-o',dest='outfile',metavar='--output',action="store",help="Output file",default=None)
parser.add_argument('-v',dest='var',metavar='--variable',action="store",help="Variable name (when using mathematical expressions)",default="x")
#parser.add_argument('-I',dest='ign',metavar='--ignore',action="store",help="Regex pattern of lines to ignore",default=None)


args=parser.parse_args().__dict__

colR=str(args["col"])
try:
    ind=colR.index("-")
    cS=int(colR[0:ind])
    cE=int(colR[ind+1:])
    fL=cE-cS+1
    fmt="<"+str(fL)+"s"
    colN=None
except ValueError:
    colN=int(colR)
fl=open(str(args["infile"]),"r")
delim=args["delim"]
expr=args["math"]
mult=float(args["mult"])
outf=None
mfirst=args["mfirst"]


if args["outfile"]: outf=open(str(args["outfile"]),"w")

for l in fl:
    l=l[0:len(l)-1] #Remove the "\n" from the end
    if colN:
        frags=str.split(l,delim)
        no=float(frags[colN-1])
    else:
        preS=l[0:cS-1]
        postS=l[cE:]
        no=float(l[cS-1:cE])

    if mfirst=="1": no*=mult
    if expr: no=float(eval(expr.replace(str(args["var"]),str(no))))
    if mfirst!="1": no*=mult

    if colN:
        frags[colN-1]=str(no)
        res=delim.join(frags)
    else:
        no=round(no,fL)
        midS=format(str(no),fmt)
        if len(midS)>fL: midS=midS[0:fL]
        res=preS+midS+postS
    if outf: outf.write(res+"\n")
    else: print(res)
if outf: outf.close()
