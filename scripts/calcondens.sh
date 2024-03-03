#!/bin/bash
# CALCULATES CONNECTION DENSTIY (NO OF CONNECTIONS IN EACH MOLECULE - EXCEPT HYDROGEN)

if test -z "stripall.sh"
then
	echo "This script requires the 'stripall.sh' script. Please copy/link it here"
	exit 2
fi
./stripall.sh
x=`ls edit*_noH.pdb`
for f in $x; do echo $f" "`grep ^CONECT $f|sed s/"  *"/" "/g|sed s/"CONECT "//g|grep -o " "|wc -l`; done > condens.dat
echo "Written to condens.dat"
