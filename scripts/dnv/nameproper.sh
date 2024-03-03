#!/bin/bash
# Rename file properly (dnv) n1_n2.pdb -> result_n1_n2.pdb

for f in `ls *.pdb`
do
	nump=`echo $f|grep -o "[0-9]*_[0-9]*.pdb"`
	if test $? -ne 0
	then
		continue
	fi
	mv $f "result_"$nump
done
