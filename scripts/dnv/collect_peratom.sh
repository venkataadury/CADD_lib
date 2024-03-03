#!/bin/bash
# Collect molecule energy to protein (dnv) and also calculate the average energy per-atom

for x in `head -1 result_prot_*.pdb|grep "==>"|cut -d" " -f2`; 
do 
	n=`sed /"H *"/d $x|grep "^ATOM"|wc -l`
	en=`head -1 $x|awk -e '{print $2}'`
	hen=`head -1 $x|awk -e '{print $3}'`
	echo $x `divi $en $n` `divi $hen $n` $n
done
echo "Written to energies.dat"
