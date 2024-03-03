#!/bin/bash

fl=$1
for f in `cat $fl`
do
	echo $f
	efn=`echo $f|sed s/'.pdb'//g`
	cd $efn
	if test -f "encal.edr"
	then
		echo "Energy is calculated. Skipping rerun"
	else
		gmx_2020 grompp -f mdps/nvt_withen.mdp -c LIG_solv_emin.gro -p LIG_GMX.top -o encal.tpr -maxwarn 3
		if test $? -ne 0
		then
			echo "GROMPP failed!"
			exit 1
		fi
		gmx_2020 mdrun -rerun simrun.trr -deffnm encal -s encal.tpr -v
		if test $? -ne 0
		then
			echo "MDrun failed"
			exit 1
		fi
	fi
	echo "0"|gmx_2020 energy -f encal.edr 2> errlog
	par="0"
	par=`cat errlog|grep "....LJ-SR:[a-zA-Z]*-non-[a-zA-Z]*" -o|cut -c1-2`" $par"
	par=`cat errlog|grep "....Coul-SR:[a-zA-Z]*-non-[a-zA-Z]*" -o|cut -c1-2`" $par"
	echo $par
	echo $par|gmx_2020 energy -f encal.edr > inten
	echo "$f "`tail -2 inten|awk -e '{print $2}'` >> ../interaction_energy.spc
	cd ..
done
