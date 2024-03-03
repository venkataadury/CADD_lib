#!/bin/bash

files=$1
ligres=$2
ligsuf=$3
protsuf=$4
if test -z "$ligsuf"
then
	ligsuf="_lig"
fi
if test -z "$protsuf"
then
	protsuf="_prot"
fi
for f in $files
do
	echo $f
	fnam=`echo $f|sed s/".[A-Za-z0-9]*$"//`
	ext=`echo $f|grep -o "\..*$"`
	grep -e "$ligres" $f > $fnam"$ligsuf"$ext
	grep -e "$ligres" -v $f > $fnam"$protsuf"$ext
done
