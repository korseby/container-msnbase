#!/bin/bash
csvfile=$1   		# single csv file containing all scored candidates
outputfile=$2		# output file (pdf)
# check input
if [ -z ${csvfile+x} ] || [ -z ${outputfile+x} ]
then
	echo "csvfile: ${csvfile}"
	echo "outputfile: ${outputfile}"
	echo "'csvfile string' and 'outputfile' missing"
	echo "try: bash visualize-metfrag-multiple.sh 'csvfile' 'outputfile'"
	exit 1
fi
if [ "$(wc -l $csvfile | cut -d" " -f1)" -eq "0" ]; then 
   echo "No results found"
   exit 0
fi
if [ "$(head -n 1 $csvfile | grep -c "^parentRT,parentMZ")" -eq "0" ]; then
   echo "CSV file has no proper format"
   exit 1
fi
# add header to output file
head -n 1 $csvfile > $outputfile
# retrieve best results per m/z rt pair
scoreColumnIndex=$(head -n 1 $csvfile | tr ',' '\n' | grep -n "^Score$" | cut -d":" -f1)
OLDIFS=$IFS
IFS=$'\n'
# select top1 candidate for each parent RT and m/z
for i in $(cut -d"," -f1,2 $csvfile | tail -n +2 | sort | uniq); do
   paste <(grep "^$i" $csvfile | cut -d"," -f$scoreColumnIndex) <(grep "^$i" $csvfile) | sort -rn | cut -f2 | head -n 1 >> $outputfile
done
IFS=$OLDIFS
fileindex=0
tr -d $'\r' < $outputfile > ${outputfile}_tmp
mv -f ${outputfile}_tmp $outputfile
