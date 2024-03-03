#!/bin/bash

#DEPENDS: prepgmx.sh,doMD.sh

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
for f in $fl
do
	echo "Started with $f"
	fnam=`echo $f|sed s/'.pdb'//g`
	mkdir -p $curfol/$fnam
	ln -s $floc/$f $curfol/$fnam
	ln -s `pwd`/mdps $curfol/$fnam
	if test -f  "`pwd`"/prepis/$f"_tautomer.prepi"
	then
		ln -s `pwd`/prepis/$f"_tautomer.prepi" $curfol/$fnam
		ln -s `pwd`/frcmods/$f"_tautomer.frcmod" $curfol/$fnam
	fi
	ln -s `pwd`/acpype.py $curfol/$fnam
	ln -s `pwd`/prepgmx.sh $curfol/$fnam
	cd $curfol/$fnam
	if test -f "nvt.tpr"
	then
		echo "In $f: Minimzation complete. Moving to MD"
	else
		./prepgmx.sh $f -vm #NOTE: Change this if you DO NOT want vacuum minimization (remove '-vm' from options)
		echo "Prepared simulation"
	fi
	$curfol/doMD.sh $1 $f -r
	cd $rf
done
