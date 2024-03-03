#!/bin/bash
# Write the size of each molecule (in heavy atoms) to a file

x=`ls result_prot_*.pdb`
for f in $x;  do n=`grep ^[HA][ET][TO][AM] $f|sed /"H *$"/d|wc -l`; echo $f" "$n; done > sizes.dat
echo "Written to sizes.dat"
