#!/bin/bash

for f in `ls *.pdb`
do
	nump=`echo $f|grep -o "[0-9]*_[0-9]*.pdb"`
	if test $? -ne 0
	then
		continue
	fi
	mv $f "result_"$nump
done
