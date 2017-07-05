#!/bin/bash

# install required packages
apt-get update -y && apt-get install -y --no-install-recommends wget ca-certificates

# retrieve test case data
mkdir /tmp/testfiles
wget -O /tmp/testfiles/PE_PG_AutoMS_2-B,3_03_18583.mzML https://github.com/phnmnl/container-MSnbase/raw/develop/testfiles/PE_PG_AutoMS_2-B%2C3_03_18583.mzML

# perform test case
/usr/local/bin/readMS2MSnBase.r input=/tmp/testfiles/PE_PG_AutoMS_2-B,3_03_18583.mzML output=/tmp/msnbase-read-msms-data.rdata

# check whether files were created
if [ ! -f /tmp/msnbase-read-msms-data.rdata ]; then exit 1; fi

# remove files
rm /tmp/msnbase-read-msms-data.rdata
