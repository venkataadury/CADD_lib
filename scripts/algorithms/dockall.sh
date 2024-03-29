#!/bin/bash
#
# Dock all molecules in a folder to a given target (PDBQT format) using AutoDock VINA


prot=$1
cx=$2
cy=$3
cz=$4
sx=$5
sy=$6
sz=$7
of=$8
hyd=$9


x=`ls *.pdb`
for f in $x
do
	if test -f $f".pdbqt"
	then
		echo "Skipping already converted file: $f"
	else
		if test -n "$hyd"
		then
			obabel -i pdb $f -o pdbqt -xh -O $f".pdbqt"  -r
		else
			obabel -i pdb $f -o pdbqt -O $f".pdbqt" -r
		fi
	fi
	if test -f $f"_out.pdbqt"
	then
		echo "Skipping docking: Result already exists"
	else
		vina --receptor $prot --ligand $f".pdbqt" --center_x $cx --center_y $cy --center_z $cz --size_x $sx --size_y $sy --size_z $sz
	fi
done
for f in $x
do
	echo $f".pdbqt"
	grep "REMARK VINA RESULT" $f"_out.pdbqt"
	if test -n "$of"
	then
		echo $f".pdbqt" >> $of
		grep "REMARK VINA RESULT" $f"_out.pdbqt" >> $of
	fi
done
grep "^REMARK VINA" *_out.pdbqt -m1 |awk -F":" '{print $1" "$3}'|awk '{print $1" "$2}' > docked_energies.dat
echo "Written energies to docked_energies.dat"
