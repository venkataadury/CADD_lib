#!/bin/bash

files=`ls *.pdb`
for f in $files
do
	./rewrite $f
done
