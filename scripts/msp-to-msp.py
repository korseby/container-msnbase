#!/usr/bin/env python3

# Load modules
import errno
import sys
import os
import argparse
import glob
from pathlib import Path
import re

# Parse arguments
parser = argparse.ArgumentParser(description='Merge several msp files into one msp and perform validation on the resulting msp file.')
parser.add_argument('-v', '--version', action='version', version='msp-to-msp Version 0.1',
                   help='show version')
parser.add_argument('-i', '--input_dir', metavar='dir', dest='input_dir', required=True,
                   help='directory containing the msp files')
parser.add_argument('-o', '--output_file', metavar='file', dest='output_file', required=True,
                   help='output msp file')
parser.add_argument('-a', '--fix-alignment-id', dest='alignment_id', action='store_true', required=False,
                   help='Create a unique AlignmentID')
parser.add_argument('-m', '--sec-to-min', dest='sec_to_min', action='store_true', required=False,
                   help='Convert rt seconds to minutes')
parser.add_argument('-s', '--min-to-sec', dest='min_to_sec', action='store_true', required=False,
                   help='Convert rt minutes to seconds')

args = parser.parse_args()

# Input files
input_dir = args.input_dir

input_files = sorted( [f for f in glob.glob(''.join([input_dir, '/', '*']))] )

# Output file
output_file = args.output_file

# Fix AlignmentID?
alignment_id = args.alignment_id

# Convert seconds to minutes?
sec_to_min = args.sec_to_min

# Convert minutes to seconds?
min_to_sec = args.min_to_sec



# Merge MSP files
msp_lib = []
for msp_filename in input_files:
	msp_file = open(msp_filename, "r")
	for line in msp_file:
		line = line.rstrip('\r\n')
		msp_lib.append(line)
	msp_file.close()



# Fix Alignment-ID
if (alignment_id == True):
	count = 1
	for lnum, line in enumerate(msp_lib, start=0):
		if ("alignmentid" in line.lower()):
			line = str('AlignmentID: ' + str(count))
			count = count+1
		msp_lib[lnum] = line



# Convert seconds to minutes
if (sec_to_min == True):
	for lnum, line in enumerate(msp_lib, start=0):
		if ("RETENTIONTIME" in line.upper()):
			rt = float(re.sub(".*\: ", "", line))
			rt = rt / 60
			line = str('RETENTIONTIME: ' + str(rt))
		msp_lib[lnum] = line



# Convert minutes to seconds
if (min_to_sec == True):
	for lnum, line in enumerate(msp_lib, start=0):
		if ("RETENTIONTIME" in line.upper()):
			rt = float(re.sub(".*\: ", "", line))
			rt = rt * 60
			line = str('RETENTIONTIME: ' + str(rt))
		msp_lib[lnum] = line



# Save output MSP
msp_file = open(output_file, "w")
for line in msp_lib:
	msp_file.write(line + '\n')
msp_file.close()



