#!/bin/bash

x=$1
y=$2

if test -z "$x"
then
	exit 0
fi
ss=$x
if test -n "$y"
then
	ss=$y
fi
searchres=`googler --noprompt "$ss"|grep -E "drugbank.ca|zinc|wikipedia"|head -1`
if test -z "$searchres"
then
	echo "Drug: '$x' not found in Drugbank or Wikipedia. Maybe you should search manually including these websites. Also try ZINC"
	exit 1
fi
searchres=`echo $searchres|sed s/"^ *"//g`
echo "Found $x as/at $searchres"
wget "$searchres" -O DRUG -nv
echo $searchres|grep "drugbank.ca" 1> /dev/null 2> /dev/null
if test $? -eq 0
then
	#Result was found at drugbank
	url="http://drugbank.ca/"`grep "SMILES" DRUG|grep -E "href=\"[A-Za-z0-9_\/]*.smiles" -o|sed s/"href=\""//g`
	echo $url
	wget -q "$url" -O $x".smi"
	exit 0
fi
echo $searchres|grep "wikipedia.org" 1> /dev/null 2> /dev/null
if test $? -eq 0
then
	#Result was found at Wikipedia
	ln=`grep SMILES -i DRUG -A1|tail -1|grep -E "[/\\A-Za-z0-9+--@=\(\)\[]*][/\\A-Za-z0-9@=\(\)]*" -o`
	if test -z "$ln"
	then
		 ln=`grep SMILES -i DRUG -A1|tail -1|grep -E "[/\\A-Za-z0-9+--@=\(\)\[]*<" -o|sed s/"<"//g|sed /"^$"/d`
	fi
	echo $ln|sed s/" "//g > $x".smi"
	exit 0
fi
echo "Returned empty"
exit 2
