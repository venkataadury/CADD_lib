#!/bin/bash

i=$1 #Start number
protpath=$3 #Absolute path to protein receptor file
cx=$4
cy=$5
cz=$6
sx=$7
sy=$8
sz=$9
while test $i -le $2 #End number
do
	cd size"$i"
	ln -s $protpath .
	ln -s ~/CADD/dockall.sh .
	./dockall.sh $protpath $cx $cy $cz $sx $sy $sz
	grep "^REMARK VINA" -m1 *_out.pdbqt|awk -F":" -e '{print $1" "$3}'|awk -e '{print $1" "$2}' > docked_energies.dat
	cd ..
	echo "Processed: $i"
	i=`expr $i + 1`
done
