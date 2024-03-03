#!/bin/bash

i=$1 #Start number
#$3 is the top n
#$4 is the protein (folder) name
if test -z "$4"
then
	prot="abhd12"
else
	prot=$4
fi
while test $i -le $2 #End number
do
	x=`sort -k2 -n ../..//$prot/results/size$i/energies.dat|head -$3|cut -d" " -f1`
	for f in $x
	do
		cp ../../$prot/results/size$i/$f .
	done
	echo "Processed: $i"
	i=`expr $i + 1`
done
