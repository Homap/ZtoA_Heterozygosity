#!/usr/bin/python
import sys
from datetime import datetime

startTime = datetime.now()

inFile = sys.argv[1]
outFile = inFile + ".fa"

print ">> Opening Fastq file..."

try:
	fastq = open(inFile, 'r')
	fasta = open(outFile, 'w')
except IOError:
	print "Failed to open" + inFile
	exit(1)

canPrintLines = False # Don't want to print anything until we get to an @
for line in fastq:
	line = line.rstrip("\n")
   	if line.startswith("@") and len(line) < 30: # 30 for all species, 40 for Melopsittacus undulatus
		canPrintLines = True # We have found an @ so we can start printing lines
		line = line.replace("@", ">")
   	elif line.startswith("+"):
		canPrintLines = False # We have found a + so we don't want to print anymore
	
	if canPrintLines: # or if canPrintLines is True:
		fasta.write(line + "\n")

print datetime.now() - startTime
		
fastq.close()
fasta.close()
