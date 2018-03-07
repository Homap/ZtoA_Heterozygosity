#!/usr/bin/python
# Created by: Homa Papoli
# Description: Converts multiline FASTAs to single line FASTAs
#
# Usage: MLtoSLFasta.py <sequences.fa>
# Example: MLtoSLFasta.py mySeqs.fa

#------------------------------------------------------------------------------------------------
#================================================================================================
#Imports
import sys
import re
import gzip
from datetime import datetime

#================================================================================================
# Store the current time:
startTime = datetime.now()
# Checks if proper number of arguments are given and gives instructions for proper use
def argsCheck(numArgs):
	if len(sys.argv) < numArgs or len(sys.argv) > numArgs:
		print "Converts multiline FASTAs to single line FASTAs"
		print "By Homa Papoli---The structure of the code has been copied from the script called fastamultitoone.py of Lee Bergstrand\n"
		print "Usage: " + sys.argv[0] + "<sequences.fasta>"
		print "Examples: " + sys.argv[0] + " mySeqs.fasta"
		exit(1) # Aborts program. (exit(1) indicated that an error occured)
#================================================================================================
# Main program code:
# House keeping...
argsCheck(2) # Checks if the number of arguments is correct.

# Input and output
inFile = sys.argv[1]
outFile = inFile + ".SL"

print ">> Opening FASTA file..."

try:
	f = open(inFile, 'r')
	o = open(outFile, 'w')
except IOError:
	print "Failed to open" + inFile
	exit(1)

print ">> Converting FASTA file from multiline to single line and writing to file"

currentline = ""
for line in f:
	if line.startswith('>'):
		line = line.rstrip('\n')
		if currentline != "": o.write(currentline + '\n')
		o.write(line + '\n')
		currentline = ""
	else:
		line = line.rstrip('\n')
		currentline = currentline + line
		
o.write(currentline + '\n')
f.close()

print ">> Done"
print datetime.now() - startTime
