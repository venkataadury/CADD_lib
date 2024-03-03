#!/bin/bash

name=$1
wget "https://www.drugbank.ca/unearth/q?utf8=✓&searcher=drugs&query=""$name" -O query.html -a dnld.log
ids=`grep -E -o -m1 "DB[0-9][0-9]*" query.html|uniq`
idc=`echo $ids|wc -w`
echo "$name~$ids"
if test $idc -eq 1
then
	a=
else 
	echo "Ambiguous ID for $name"
fi
