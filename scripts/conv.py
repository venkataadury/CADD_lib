#!/usr/bin/python3

import sys
import os

fl=open(sys.argv[1],"r")
for l in fl:
    l=l[0:len(l)-1] #Remove the '\n' from the end
    cols=str.split(l,' ')
    smstr=cols[0]
    flnm=cols[1]
    print("echo '"+smstr+"'| obabel -ismi -o"+sys.argv[2]+" -O "+flnm+"."+sys.argv[2])
    os.system("echo '"+smstr+"'| obabel --gen3d -ismi -o"+sys.argv[2]+" -O "+flnm+"."+sys.argv[2])
