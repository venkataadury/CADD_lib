#!/bin/bash

# DeNovo push files. A quick and dirty way to increment all file name IDs (numerical) by a constant
shiftval=$1
templates=$2
field=$3
if test -z "$templates"
then
	template="result_" #.pdb, .mol2
fi
if test -z "$field"
then
	field=0
fi

echo "Starting move of file. DO NOT INTERRUPT THIS PROCESS!!!"
read -p "Press enter to continue"
randfol=`uuidgen -r`
mkdir -p "$randfol"
cd "$randfol"
echo "Example:"
find .. -maxdepth 1 -name "$templates*.pdb"|head -1|sed s/"..\/"//g|xargs -I"%" sh -c "echo ../% && numarith % + $shiftval $field"|xargs -n 2 -p mv
find .. -maxdepth 1 -name "$templates*.pdb"|sed s/"..\/"//g|xargs -I"%" sh -c "echo ../% && numarith % + $shiftval $field"|xargs -n 2 mv
cd ..
find "$randfol" -name "*.pdb"|xargs -I"%" mv % .
rmdir "$randfol" || echo "WARN: Couldn't delete temporary directory" >&2
echo "Completed move"
