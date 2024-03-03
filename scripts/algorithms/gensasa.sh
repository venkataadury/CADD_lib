#!/bin/bash
#
# Generate SASA values for a ligand by simulation

GENNAME=LIG
echo "Using general name $GENNAME"
TOPFILE=$GENNAME"_GMX.top"
GROFILE=$GENNAME"_GMX.gro"
outf=$1
if test -z "$outf"
then
	outf="SASA_results.spc"
fi
echo "Remember to run in amber mode"
for f in `ls *.pdb`
do
	echo "Started $f"

	if test -f "$f""_tautomer.prepi"
	then
		echo "File $f was processed by antechamber already. Skipping."
		grep "^$f" $outf
		if test $? -eq 0
		then
			echo "Result found. Continuing with next molecule"
			continue
		fi
	else
		antechamber -fi pdb -fo prepi -i "$f"  -o $f"_tautomer.prepi" -c bcc -nc 0 -rn $GENNAME
		if test $? -ne 0
		then
			echo "antechamber stopped abruptly."
			echo "Exiting..."
			exit 1
		fi
	fi
	if test -f $f"_topol.tpr"
	then
		echo
	else
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
	gmx editconf -f $GROFILE -box 3 -bt cubic -d 1 -center 0 -o LIG.gro
	if test $? -ne 0
	then
		echo  "Editconf failed"
		exit 2
	fi
	gmx solvate -cp LIG.gro -cs spc216 -p $TOPFILE -o LIG_solv.gro
	if test $? -ne 0
	then
		echo "Solvate failed"
		exit 2
	fi
	head -5 $TOPFILE > temp.top
	echo "#include \"amber99sb.ff/ffnonbonded.itp\"" >> temp.top
	tail +6 $TOPFILE|head -n -7 >> temp.top
	echo "#include \"amber99sb.ff/spc.itp\"" >> temp.top
	tail -7 $TOPFILE >> temp.top
	cp $TOPFILE backup.top
	mv temp.top $TOPFILE
	gmx grompp -f mdps/test.mdp -c LIG_solv.gro -p LIG_GMX.top -o $f"_topol.tpr"
	if test $? -ne 0
	then
		echo "WARN: grompp failed"
	fi
	test -f $f"_topol.tpr"
	if test $? -ne 0
	then
		echo "$f""_topol.tpr not found when processing '$f'"
		exit 3
	fi
	fi
	echo "2"|gmx sasa -s $f"_topol.tpr" -o test
	echo "$f `tail -1 test.xvg`" >> $outf
	rm -f test.xvg LIG.gro LIG_solv.gro $TOPFILE $GROFILE
done
