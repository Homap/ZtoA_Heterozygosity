#!/usr/bin/python
# Created by: Homa Papoli
# Description: Compares the contents of the sequences of two FASTA files
#
# Usage: compare_scaffolds.py <seq1.fa> <seq2.fa>

#-----------------------------------------------------------------------
# Imports
import sys
from FastaModule import *
#from datetime import datetime

#-----------------------------------------------------------------------
# Store the current time
#startTime = datetime.now()
## Checks if proper number of arguments are given and gives instructions for proper use
#def argsCheck(numArgs):
#	if len(sys.argv) < numArgs or len(sys.argv) > numArgs:
#		print "Compare the contents of the sequences of two FASTA files"
#		print "Usage: " + sys.argv[0] + "<seq1.fa> + <seq2.fa>"
#		exit(1)
##-----------------------------------------------------------------------
# Main program code:
# House keeping...
#argsCheck(4) # Check if the number of arguments is correct.
		
# Input and output
try:
	fasta1 = open(sys.argv[1], 'r')
	fasta2 = open(sys.argv[2], 'r')
	output = open(sys.argv[3], 'w')
except IOError as ex:
	print "Sorry, failed to open: " + ex.strerror
	exit(1)
# Read Fasta files into dictionaries
print 'Reading fasta1 into dictionary'
fastaDict1 = FaAsDict(fasta1)
print 'Reading fasta2 into dictionary'
fastaDict2 = FaAsDict(fasta2)
	
#{keys.split()[0]: value for keys, value in fastaDict2.items()}
# In the fasta genome assembly file, the header of some sequences consist of the name and of a number. I change the header by getting only the first part of the name.
for keys in fastaDict2.keys():
	fastaDict2[keys.split()[0]] = fastaDict2.pop(keys)	
	
for keys in fastaDict1.keys():
	fastaDict1[keys] = fastaDict1[keys].upper()
#	print fastaDict1
		
for keys in fastaDict2.keys():
	fastaDict2[keys] = fastaDict2[keys].upper()
#	print fastaDict2		

IUPAC = ['R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V', 'N']

#for keys in fastaDict1.keys():
#	output.write(keys + '\n')
#	if keys in fastaDict2.keys():
#		for base in range(0, len(fastaDict1[keys])):
#			if fastaDict1[keys][base] in IUPAC:#== 'N':
#				pass
#			elif fastaDict1[keys][base] != 'N' and fastaDict1[keys][base] == fastaDict2[keys][base]:
#				pass
#			else:
#				output.write(fastaDict1[keys][base] + '\t' + fastaDict2[keys][base] + '\t' + str(base) + '\n')	
seqs = []				
for keys in fastaDict1.keys():
	print keys
	if keys in fastaDict2.keys():
#		for base in range(0, len(fastaDict1[keys])):
#			if fastaDict1[keys][base] in IUPAC:#== 'N':
#				pass
		if fastaDict1[keys][len(fastaDict1[keys])-20 : len(fastaDict1[keys])] == fastaDict2[keys][len(fastaDict1[keys])-20 : len(fastaDict1[keys])]:
			pass
		else:
			seqs = [fastaDict1[keys][len(fastaDict1[keys])-20 : len(fastaDict1[keys])], fastaDict2[keys][len(fastaDict1[keys])-20 : len(fastaDict1[keys])]]# keys + '\t' + 'Last 20 bases do not match!'	
		print seqs
#		for base in range(0, len(seqs[0])):
#			if seqs[0][base] != seqs[1][base]:
#				print seqs[0][base] + '\t' + seqs[1][base] + '\t' + str(base)
	
#	print seqs			
					

fasta1.close()
fasta2.close()








