#!/bin/sh

#  TransformModelCheckingOutput.sh
#
# $1 Path to Prism-output file

awk '/^[0123456789]+\:/' $1 > tmp.txt
sed -i "" -E 's/[0-9]+\://g' tmp.txt
sed -i "" 's/Infinity/NA/g' tmp.txt
sed -i "" 's/(//g' tmp.txt
sed -i "" 's/)//g' tmp.txt
sed -i "" 's/=/    /g' tmp.txt
sed -i "" 's/,/    /g' tmp.txt
mv tmp.txt $1
