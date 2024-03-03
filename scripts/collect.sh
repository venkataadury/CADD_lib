#!/bin/bash

of=$1
if test -z $of
then
	of="energies.dat"
fi
for x in `head -1 result_prot_*.pdb|grep "==>"|cut -d" " -f2`; do echo $x `head -1 $x|cut -d" " -f2`; done > $of
n=`cat $of|wc -l`
if test $n -gt 0
then
	echo "Written to $of"
else
	echo "No files found"
	rm -rf $of
fi
