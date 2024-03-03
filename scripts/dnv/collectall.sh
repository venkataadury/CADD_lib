#!/bin/bash

i=$1 #Start number
while test $i -le $2 #End number
do
	cp ../size$i/result*.pdb .
	echo "Processed: $i"
	i=`expr $i + 1`
done
