#!/bin/bash
csvfiles=$1   		# comma separated list of MetFrag csv files
outputfile=$2		# output file (pdf)
parameterfiles=$3	# file containing commandline parameters used to generate files
# check input
if [ -z ${csvfiles+x} ] || [ -z ${outputfile+x} ]
then
	echo "csvfiles: ${csvfiles}"
	echo "outputfile: ${outputfile}"
	echo "'csvfiles string' and 'outputfile' missing"
	echo "try: bash visualize-metfrag.sh 'csvfiles' 'outputfile'"
	exit 1
fi
# add header to output file
for i in $(echo $csvfiles | tr ',' '\n')
do
	head -n 1 $i
done | head -n 1 > $outputfile
fileindex=0
if [ ! -z ${parameterfiles+x} ]	&& [ "${parameterfiles}" != "" ]						# add url strings to result files if set
then 
	IFS=',' read -a paramArray <<< "$parameterfiles"
	sed -i "s/\(.*\)/\1,SampleName,MetFragWebURL/" $outputfile		# add additional columns
	for i in $(echo $csvfiles | tr ',' '\n') 
	do
		line=$(cat ${paramArray[${fileindex}]})
		numLines=$(wc -l $i | cut -d' ' -f1)
		if [ "$numLines" -gt "1" ]			# only add parameters if csv file has entries
		then
			# convert the parameters 
			sampleName=$(echo $line | sed "s/.*SampleName\s*=\s*//" | sed "s/\s.*//" | sed "s/.*\/\(.*\)\.txt/\1/")		# get SampleName
			url="http://msbi.ipb-halle.de/MetFragBeta/landing.xhtml"
			# fetch parameters for MetFragWeb URL
			FragmentPeakMatchAbsoluteMassDeviation=$(echo $line | tr ' ' '\n' | grep "^FragmentPeakMatchAbsoluteMassDeviation=" | cut -d"=" -f2)
			FragmentPeakMatchRelativeMassDeviation=$(echo $line | tr ' ' '\n' | grep "^FragmentPeakMatchRelativeMassDeviation" | cut -d"=" -f2)
			DatabaseSearchRelativeMassDeviation=$(echo $line | tr ' ' '\n' | grep "^DatabaseSearchRelativeMassDeviation=" | cut -d"=" -f2)
			IonizedPrecursorMass=$(echo $line | tr ' ' '\n' | grep "^PrecursorCompoundIDs=" | cut -d"=" -f2)
			NeutralPrecursorMass=$(echo $line | tr ' ' '\n' | grep "^NeutralPrecursorMass=" | cut -d"=" -f2)
			NeutralPrecursorMolecularFormula=$(echo $line | tr ' ' '\n' | grep "^NeutralPrecursorMolecularFormula=" | cut -d"=" -f2)
			PrecursorIonMode=$(echo $line | tr ' ' '\n' | grep "^PrecursorIonMode=" | cut -d"=" -f2)
			PeakList=$(echo $line | tr ' ' '\n' | grep "^PeakListString=" | cut -d"=" -f2)
			MetFragDatabaseType=$(echo $line | tr ' ' '\n' | grep "^MetFragDatabaseType=" | cut -d"=" -f2)
			# start adding parameters to url
			if [ "$MetFragDatabaseType" == "MetChem" ]; then MetFragDatabaseType=KEGG; fi
			url=${url}?MetFragDatabaseType=$MetFragDatabaseType
			if [ "$FragmentPeakMatchAbsoluteMassDeviation" != "" ]; then url="${url}&FragmentPeakMatchAbsoluteMassDeviation=$FragmentPeakMatchAbsoluteMassDeviation"; fi
			if [ "$FragmentPeakMatchRelativeMassDeviation" != "" ]; then url="${url}&FragmentPeakMatchRelativeMassDeviation=$FragmentPeakMatchRelativeMassDeviation"; fi
			if [ "$DatabaseSearchRelativeMassDeviation" != "" ]; then url="${url}&DatabaseSearchRelativeMassDeviation=$DatabaseSearchRelativeMassDeviation"; fi
			if [ "$IonizedPrecursorMass" != "" ]; then url="${url}&IonizedPrecursorMass=$IonizedPrecursorMass"; fi
			if [ "$NeutralPrecursorMass" != "" ]; then url="${url}&NeutralPrecursorMass=$NeutralPrecursorMass"; fi
			if [ "$NeutralPrecursorMolecularFormula" != "" ]; then url="${url}&NeutralPrecursorMolecularFormula=$NeutralPrecursorMolecularFormula"; fi
			if [ "$PrecursorIonMode" != "" ]; then url="${url}&PrecursorIonMode=$PrecursorIonMode"; fi
			if [ "$PeakList" != "" ]; then url="${url}&PeakList=$PeakList"; fi
			echo "$(head -n 2 $i | tail -n 1),\"${sampleName}\",\"$url\"" >> $outputfile			
			# to check
		fi
		fileindex=$((fileindex+1))
	done
else									       # just get first result entry
	for i in $(echo $csvfiles | tr ',' '\n')
        do
		numLines=$(wc -l $i | cut -d' ' -f1)
               	if [ "$numLines" -ne "0" ]                   # if file is not empty
               	then
               	        # convert the parameters 
                       	head -n 2 $i | tail -n 1 >> $outputfile                    
               	fi
	done
fi
tr -d $'\r' < $outputfile > ${outputfile}_tmp
mv -f ${outputfile}_tmp $outputfile
