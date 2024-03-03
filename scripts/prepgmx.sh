#!/bin/bash

GENNAME=LIG
echo "Using general name $GENNAME"
TOPFILE=$GENNAME"_GMX.top"
GROFILE=$GENNAME"_GMX.gro"
f=$1
vm=$2
antec=$3
if [ "$vm" == "--ac" -o "$vm" == "-ac" ]; then
	vm=$antec
	antec="--ac"
fi
if test -z "$f"
then
	echo "No input file provided"
	exit 127
fi
echo "Remember to run in amber mode"
	echo "Started $f"

	if test -f "$f""_tautomer.prepi"
	then
		echo "File $f was processed by antechamber already. Skipping."
	else
		sed s/"FA"/"F1"/g $f > temp.pdb
		mv temp.pdb $f
		antechamber -fi pdb -fo prepi -i "$f"  -o $f"_tautomer.prepi" -c bcc -nc 0 -rn $GENNAME
		if test $? -ne 0
		then
			echo "antechamber stopped abruptly."
			echo "Exiting..."
			exit 1
		fi
	fi
	parmchk2 -i $f"_tautomer.prepi" -o tautomer.frcmod -f prepi

	echo -e "loadamberparams tautomer.frcmod\nloadamberprep $f""_tautomer.prepi\nlist\ncheck $GENNAME unit\nsaveamberparm $GENNAME tautomer.prmtop tautomer.inpcrd\nquit" > leaprc.leapin
	tleap -s -f leaprc.gaff -f leaprc.leapin
	python2.7 acpype.py -p tautomer.prmtop -x tautomer.inpcrd
	if test $? -ne 0
	then
		echo "acpype failed"
		exit 4
	fi
	echo "Processed $f"
	gmx_2020 editconf -f $GROFILE -box 3.75 -bt cubic -c -d 1 -o LIG.gro
	if test $? -ne 0
	then
		echo  "Editconf failed"
		exit 2
	fi
	head -5 $TOPFILE > temp.top
	echo "#include \"amber99sb.ff/ffnonbonded.itp\"" >> temp.top
	tail +7 $TOPFILE|head -n -7 >> temp.top
	echo "#include \"amber99sb.ff/spc.itp\"" >> temp.top
	tail -7 $TOPFILE >> temp.top
	cp $TOPFILE backup.top
	mv temp.top $TOPFILE
	if test -n "$vm" #vm=='--vm'
	then
		if test -f "vacmin.gro"
		then
			echo "Vacuum minimization complete"
			cp vacmin.gro LIG.gro
		else
			echo "Preforming vacuum minimization"
			gmx_2020 grompp -f mdps/vacmin.mdp -c LIG.gro -p LIG_GMX.top -o vacmin.tpr -maxwarn 3
			if test $? -eq 0
			then
				gmx_2020 mdrun -deffnm vacmin -s vacmin.tpr -v
				if test $? -eq 0
				then
					mv LIG.gro LIG_orig.gro
					cp vacmin.gro LIG.gro
				else
					echo "WARN: Skipping vacuum minimzation result collection. MDrun failed"
				fi
			else
				echo "WARN: Skipping vacuum minimization step. grompp failed."
			fi
		fi
	fi
	gmx_2020 solvate -cp LIG.gro -cs spc216 -p $TOPFILE -o LIG_solv.gro
	if test $? -ne 0
	then
		echo "Solvate failed"
		exit 2
	fi
	if test -n "$antec" #Last parameter is '--ac'
	then
		exit 0
	fi
	if test -f "LIG_solv_emin.gro"
	then
		echo "Energy minimization (post solvation) complete. Skipping to MD"
	else
		gmx_2020 grompp -f mdps/minim.mdp -c LIG_solv.gro -p LIG_GMX.top -maxwarn 3
		if test $? -ne 0
		then
			echo "WARN: grompp failed"
		fi
		test -f topol.tpr
		if test $? -ne 0
		then
			echo "topol.tpr not found when processing '$f'"
			exit 3
		fi
		echo "Ready for GROMACS simulation"
		echo "Starting emin"
		gmx_2020 mdrun -deffnm LIG_solv_emin -s topol.tpr -v
		if test $? -eq 0
		then
			echo "Done! Output of energy minimization written to 'LIG_solv_emin.gro'"
		else
			echo "$f: Energy minimzation terminated improperly. Exiting..."
			exit 5
		fi
	fi
	echo "Preparing MD ... "
	gmx_2020 grompp -f mdps/nvt.mdp -c LIG_solv_emin.gro -p LIG_GMX.top -o nvt.tpr -maxwarn 3
	if test $? -eq 0
	then
		echo "Prepared successfully"
	else
		echo "Failed MD simulation preparation"
		exit 5
	fi
