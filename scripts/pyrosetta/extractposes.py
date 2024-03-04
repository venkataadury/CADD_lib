#!/usr/bin/env python3
from pyrosetta import init
init()

from pyrosetta.io import poses_from_silent
import sys

out_template=sys.argv[1]
out_template=str.split(out_template,".")[0]
print("Saving to PDB files with output template:",out_template)
for i,pose in enumerate(poses_from_silent("rosetta.out")):
    pose.dump_pdb(out_template+"_"+str(i)+".pdb")
    print("Pose",i+1,"completed",flush=True)
print("Done",flush=True)
