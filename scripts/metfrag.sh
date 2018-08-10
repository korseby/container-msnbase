#!/bin/bash
while [[ $# -gt 1 ]]
do
key="$1"

case $key in
    -z|--zip)
    ZIPFILE="$2"
    shift # past argument
    ;;
    -p|--parameters)
    PARAM="$2"
    shift # past argument
    ;;
    -d|--database)
    DB="$2"
    shift # past argument
    ;;
    -o|--output)
    OUTPUT="$2"
    shift # past argument
    ;;
    --default)
    DEFAULT=YES
    ;;
    *)
    ;;
esac
shift # past argument or value
done
echo FILE EXTENSION  = "${ZIPFILE}"
echo SEARCH PATH     = "${PARAM}"
echo LIBRARY PATH    = "${DB}"

mkdir -p tmpmetf
mkdir -p resultsmet
unzip "${ZIPFILE}" -d tmpmetf

counter=0
for filename in tmpmetf/*.*; do
counter=$((counter+1))
FIleOutName=$(basename ${filename%.*}.csv)


/usr/local/bin/metfrag -Xmx2048m -Xms1024m "`cat $filename` ResultsFile=resultsmet/${FIleOutName} NumberThreads=1 $PARAM " LocalDatabasePath=${DB}&

if (( $counter % 10 == 0 )); then wait; fi # Limit to 10 jobs at the same time

done
wait
zip -r -j metfragres.zip resultsmet/*.*

cp metfragres.zip ${OUTPUT}




