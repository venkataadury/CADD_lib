#!/bin/bash
# Parameterize a ligand from MD using antechamber

curfol=`pwd`
cd $1
rf=`pwd`
floc=$rf
if test -n "$2"
then
	cd $curfol
	cd $2
	floc=`pwd`
fi
fl=`ls *.pdb`
cd $rf
. $curfol/gromfin.sh
for f in $fl
do
	echo "Started with $f"
	fnam=`echo $f|sed s/'.pdb'//g`
	if test -d "$curfol/$fnam"
	then
		echo "Folder exists for '$f'. Skipping the processing step"
		continue
	fi
	mkdir -p $curfol/$fnam
	ln -s $floc/$f $curfol/$fnam
	ln -s `pwd`/mdps $curfol/$fnam
	if test -f  "`pwd`"/prepis/$f"_tautomer.prepi"
	then
		ln -s `pwd`/prepis/$f"_tautomer.prepi" $curfol/$fnam
	fi
	ln -s `pwd`/acpype.py $curfol/$fnam
	#ln -s `pwd`/prepgmx.sh $curfol/$fnam
	cd $curfol/$fnam
	preplig $f
	if test $? -ne 0
	then
		echo "Antechamber failed to parametrize '$f'. Deleting folder"
		cd ..
		rm -r $fnam
	fi
	cd $rf
done
