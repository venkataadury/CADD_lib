#!/bin/bash

if test -z "$GENNAME"
then
	export GENNAME=LIG
fi
echo "Using general name $GENNAME"
export TOPFILE=$GENNAME"_GMX.top"
export GROFILE=$GENNAME"_GMX.gro"
export MAXWARN=3

function preplig()
{
	f=$1
	rnc=$2
	if test -z "$rnc"
	then
		rnc="0"
	fi
	if test -f "$f""_tautomer.prepi"
	then
		echo "File $f was processed by antechamber already. Skipping."
	else
		sed s/"FA"/"F1"/g $f > temp.pdb
		mv temp.pdb $f
		antechamber -fi pdb -fo prepi -i "$f"  -o $f"_tautomer.prepi" -c bcc -nc $rnc -rn $GENNAME
		if test $? -ne 0
		then
			echo "antechamber stopped abruptly."
			echo "Exiting..."
			return 1
		fi
	fi
	parmchk2 -i $f"_tautomer.prepi" -o tautomer.frcmod -f prepi

	echo -e "loadamberparams tautomer.frcmod\nloadamberprep $f""_tautomer.prepi\nlist\ncheck $GENNAME unit\nsaveamberparm $GENNAME tautomer.prmtop tautomer.inpcrd\nquit" > leaprc.leapin
	tleap -s -f leaprc.gaff -f leaprc.leapin
	python2.7 acpype.py -p tautomer.prmtop -x tautomer.inpcrd
	if test $? -ne 0
	then
		echo "acpype failed"
		return 4
	fi
	echo "Processed $f using antechamber"
	if test -f $GENNAME"_GMX.top"
	then
		echo "Topology found. Antechamber was successful."
	else
		echo $GENNAME"_GMX.top (Topology) not found!!!"
		return 4
	fi
	#echo "Generating in-place."
	#obabel -ipdb -ogro -O $GENNAME"_GMX.gro" NEWPDB.PDB -xn
	#echo "Done."
}
function boxify()
{
	sz=$1
	additps=$2
	if test -z "$sz"
	then
		sz=3.75
	fi
	gmx_2020 editconf -f $GROFILE -box $sz -bt cubic -c -d 1 -o $GENNAME".gro"
	if test $? -ne 0
	then
		echo  "Editconf failed"
		return 2
	fi
	
	head -5 $TOPFILE > temp.top
	if test -z "$additps"
	then
		echo "#include \"amber99sb.ff/ffnonbonded.itp\"" >> temp.top
		tail +6 $TOPFILE|head -n -7 >> temp.top
		echo "#include \"amber99sb.ff/spc.itp\"" >> temp.top
		echo "#include \"amber99sb.ff/ions.itp\"" >> temp.top
		echo >> temp.top
		tail -7 $TOPFILE >> temp.top
		cp $TOPFILE backup.top
		mv temp.top $TOPFILE
	else
		cp $TOPFILE backup.top
		cp $TOPFILE temp.top
	fi
	#awk -e '{if(/; *O[A-Za-z0-9]* *- *P[L0-9] *- *O[A-Za-z0-9]* *$/) {print substr($0,0,31)"1.0947e+02"substr($0,42)} else {print $0}}' $TOPFILE > temp.top
	awk -e '{if(/; *O[A-Za-z0-9]* *- *P[L0-9] *- *O[A-Za-z0-9]* *$/) {print substr($0,0,31)"1.0947e+02"substr($0,42,13)"4"substr($0,56)} else {print $0}}' $TOPFILE > temp.top
	cp $TOPFILE $GENNAME"_GMX_orig.top"
	mv temp.top $TOPFILE
}
function vacmin()
{
	mdpfile=$1
	if test -z "$mdpfile"
	then
		mdpfile="vacmin.mdp"
	fi
	if test -f "vacmin.gro"
	then
		echo "Vacuum minimization complete"
		cp vacmin.gro $GENNAME".gro"
	else
		echo "Preforming vacuum minimization"
		if test -n "$2"
		then
			gmx_2020 grompp -f mdps/$mdpfile -c $GENNAME".gro" -p $TOPFILE -o vacmin.tpr -maxwarn $MAXWARN -r $GENNAME".gro"
		else
			gmx_2020 grompp -f mdps/$mdpfile -c $GENNAME".gro" -p $TOPFILE -o vacmin.tpr -maxwarn $MAXWARN
		fi
		if test $? -eq 0
		then
			gmx_2020 mdrun -deffnm vacmin -s vacmin.tpr -v
			if test $? -eq 0
			then
				mv $GENNAME".gro" $GENNAME"_orig.gro"
				cp vacmin.gro $GENNAME".gro"
			else
				echo "WARN: Skipping vacuum minimzation result collection. MDrun failed"
			fi
		else
			echo "WARN: Skipping vacuum minimization step. grompp failed."
		fi
	fi
}
function solvate()
{
	gmx_2020 solvate -cp $GENNAME.gro -cs spc216 -p $TOPFILE -o $GENNAME"_solv.gro"
	if test $? -ne 0
	then
		echo "Solvate failed"
		exit 2
	fi
}
function fullemin()
{
	mdpfile=$1
	if test -z "$mdpfile"
	then
		mdpfile="minim.mdp"
	fi
	if test -f $GENNAME"_solv_emin.gro"
	then
		echo "Energy minimization (post solvation) complete. Skipping to MD"
	else
		if test -n "$2"
		then
			gmx_2020 grompp -f mdps/$mdpfile -c $GENNAME"_solv.gro" -p $TOPFILE -maxwarn $MAXWARN -r $GENNAME"_solv.gro"
		else
			gmx_2020 grompp -f mdps/$mdpfile -c $GENNAME"_solv.gro" -p $TOPFILE -maxwarn $MAXWARN
		fi
		if test $? -ne 0
		then
			echo "WARN: grompp failed"
			return 1
		fi
		test -f topol.tpr
		if test $? -ne 0
		then
			echo "topol.tpr not found when processing."
			return 2
		fi
		echo "Ready for GROMACS simulation"
		echo "Starting emin"
		gmx_2020 mdrun -deffnm $GENNAME"_solv_emin" -s topol.tpr -v
		if test $? -eq 0
		then
			echo "Done! Output of energy minimization written to '$GENNAME""_solv_emin.gro'"
		else
			echo "$f: Energy minimzation terminated improperly. Exiting..."
			return 2
		fi
	fi
}
function prepMD()
{
	mdpfile=$1
	if test -z "$mdpfile"
	then
		mdpfile="nvt.mdp"
	fi
	echo "Preparing MD ... "
	if test -n "$2"
	then
		gmx_2020 grompp -f mdps/$mdpfile -c $GENNAME"_solv_emin.gro" -p $TOPFILE -o nvt.tpr -maxwarn $MAXWARN -r "$GENNAME""_solv_emin.gro"
	else
		gmx_2020 grompp -f mdps/$mdpfile -c $GENNAME"_solv_emin.gro" -p $TOPFILE -o nvt.tpr -maxwarn $MAXWARN
	fi
	if test $? -eq 0
	then
		echo "Prepared successfully for MD (nvt)"
		return 0
	else
		echo "Failed MD simulation preparation"
		return 1
	fi
}
function doMD()
{
	gmx_2020 mdrun -deffnm simrun -v -s nvt.tpr
	if test $? -ne 0
	then
		echo "ERR: MDrun failed!"
		return 1
	else
		return 0
	fi
}

#Wrappers
function doProperMD()
{
	doMD
	if test $? -ne 0
	then
		echo "MD failed! Will attempt further energy minimization"
		mv $GENNAME"_solv_emin.gro" $GENNAME"_solv.gro"
		fullemin "furtherminim.mdp"
		doMD
		if test $? -ne 0
		then
			echo "All attempts at MD failed!"
			return 2
		fi
	fi
	return 0
}
function completeProcess()
{
	f=$1
	vm=$2
	if test -f "simrun.gro"
	then
		echo "MD completed already... Skipping $f"
		return 127
	fi
	preplig $1
	boxify
	if test $? -ne 0
	then
		echo "Stopping after this..."
		return 4
	fi
	if test -n "$vm"
	then
		vacmin "vacmin.mdp"
	fi
	solvate
	fullemin
	if test $? -ne 0
	then
		echo "Could not complete energy minimization of system!"
		return 1
	fi
	prepMD
	doProperMD
	return $?
}

#Protein-Ligand
function mergesystems()
{
	x=$1
	y=$2
	t1=$3
	t2=$4
	z=$5
	to=$6
	if test -z "$z"
	then
		z="merged.gro"
	fi
	if test -z "$to"
	then
		tx=`echo $z|sed s/".gro"//g`
		to=$tx".top"
	fi
	pwd
	s1=`head -2 $x|tail -1`
	s2=`head -2 $y|tail -1`
	s=`expr $s1 "+" $s2`
	echo "Merged System" > $z
	echo " $s" >> $z
	head -n -1 "$x"|tail -n +3 >> $z
	tail -n +3 "$y" >> $z
	if test -z "$t2"
	then
		echo "No topology provided. Not generating merged topology"
		exit 0
	fi
	echo $t2|grep ".top$" 1> /dev/null 2> /dev/null
	if test $? -eq 0
	then
		tex=`echo $t2|sed s/".top"/".itp"/g`
		awk -e 'BEGIN {x=0; y=1} {if(/system/) x=1; if(/atomtypes/) y=0; if(x==0 && y==0) print}' $t2 > $tex
		if test $? -ne 0
		then
			return 2
		fi
		echo "Generated itp from top file"
		t2=$tex
	fi
	awk -e 'BEGIN {x=0} {if(/include/ && x==0) {x=1; print; print "#include \"'$t2'\""} else {print}}' $t1 > $to
	if test $? -ne 0
	then
		return 1
	fi
	rn=`head -5 $y|tail -1|cut -c6-10`
	echo "$rn         1" >> $to
	return $?
}
function addions()
{
	gmx_2020 grompp -f mdps/minim.mdp -c $GENNAME"_solv.gro" -p $TOPFILE -o traj.tpr -maxwarn $MAXWARN
	if test $? -ne 0
	then
		echo "Trajectory preparation failed"
		return 1
	fi
	pi=$1
	ni=$2
	if test -z "$pi"
	then
		pi="NA"
	fi
	if test -z "$ni"
	then
		ni="CL"
	fi
	echo "15"|gmx_2020 genion -pname $pi -nname $ni -neutral -s traj.tpr -p $TOPFILE -o ions.gro
	if test $? -ne 0
	then
		echo "Adding ions failed"
		return 2
	fi
	mv ions.gro $GENNAME"_solv.gro"
}
function fullProteinLigand()
{
	if test -f "simrun.gro"
	then
		echo "MD completed. Skipping"
		return 127
	fi
	echo $1 $2 $3 $4
	mergesystems $1 $2 $3 $4 $GROFILE $TOPFILE
	if test $? -ne 0
	then
		echo "Merging protein with ligand failed"
		return 8
	fi
	boxify 3.5 -d
	if test -n "$5"
	then
		vacmin "vacmin_restr.mdp" -r
	fi
	solvate
	addions
	if test $? -ne 0
	then
		echo "Failed in genion (for neutralization)"
		return 1
	fi
	fullemin 
	if test $? -ne 0
	then
		echo "Could not complete energy minimization of system!"
		return 1
	fi
	prepMD "nvt_big.mdp" -r
	doProperMD
	return $?
}
#fullProteinLigand $1 $2 $3 $4
#completeProcess $1 $2
