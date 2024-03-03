#!/bin/bash

i=$1 #Start number
#$3 is the top n
while test $i -le $2 #End number
do
	x=`sort -k2 -n ../size$i/energies.dat|head -$3|cut -d" " -f1`
	for f in $x
	do
		cp ../size$i/$f .
	done
	echo "Processed: $i"
	i=`expr $i + 1`
done
