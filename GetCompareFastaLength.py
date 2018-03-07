#!/usr/bin/python

# Imports
import sys
from FastaModule import *
from datetime import datetime

#---------------------------------------------------------------------------
# Store the current time:
startTime = datetime.now()
#---------------------------------------------------------------------------
# Check if the correct number of arguments have been given
def argsCheck(NumofArg):
	if len(sys.argv) < NumofArg or len(sys.argv) > NumofArg:
		print "Takes two fasta files and compare the lengths"
		print "Written by Homa Papoli"
		print "Usage: " + sys.argv[0] + "<seq1.fa> + <seq2.fa> + <out.csv>"
		exit(1)
#---------------------------------------------------------------------------
# Main program code:
# House kepping...
argsCheck(4) # Checks if the number of arguments is correct.

# Input and output
fasta1 = sys.argv[1]
fasta2 = sys.argv[2]
outFile = sys.argv[3]
# First fasta ==========================================================
print ">> Opening the first FASTA file..."

try:
	fa1 = open(fasta1, "r")
	fa2 = open(fasta2, "r")
	out = open(outFile, "w")
except IOError as ex:
	print "Sorry, failed to open: " + ex.strerror
	exit(1)

# Defining getfalen function to get the length of a fasta sequence----------
#---------------------------------------------------------------------------
print 'Reading fasta1 into dictionary'
fastaDict1 = FaAsDict(fa1)
print 'Reading fasta2 into dictionary'
fastaDict2 = FaAsDict(fa2)
	
def getfalen(seqs):
	length = {}
	for key in seqs.keys():
		length[key] = len(seqs[key])
	return length

# Calling the function for the sequences------------------------------------
#---------------------------------------------------------------------------
print ">> Get the length for fasta1"		
length1 = getfalen(fastaDict1)
#for key, value in length1.items():
#	print key + "\t" + str(value)
	
print ">> Get the length for fasta2"	
length2 = getfalen(fastaDict2)
#for key, value in length2.items():
#	print key + "\t" + str(value)
	
print ">> Compare the lengths of sequences between two fastas"
for key in length1.keys():
	if key in length2.keys() and length1[key] != length2[key]:
		out.write(key + "\t" + str(length1[key]) + "\t" + str(length2[key]) + '\n')
		
fa1.close()
fa2.close()

print ">> Done"
print datetime.now() - startTime
	
