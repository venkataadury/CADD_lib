#!/bin/bash
# Copy CONECT records from a PDBs in a folder of old-files to PDBs in a folder of new files (which lack these records)
# The script only pastes the CONECT records. It does not check them

x=$1 #Folder of old files
y=$2 #Folder of new files

ofiles=`ls $x/*.pdb`
for f in $ofiles
do
	fnam=`basename $f|cut -d"." -f1`
	grep ^CONECT $f >> $y/$fnam"_rewritten.pdb"
	echo Processed $f
done
