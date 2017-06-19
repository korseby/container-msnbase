#!/bin/bash

# perform test case
/usr/local/bin/readMS2MSnBase.r input=/testfiles/PE_PG_AutoMS_2-B,3_03_18583.mzML output=/tmp/msnbase-read-msms-data.rdata

# check whether files were created
if [ ! -f /tmp/msnbase-read-msms-data.rdata ]; then exit 1; fi

# remove files
rm /tmp/msnbase-read-msms-data.rdata
