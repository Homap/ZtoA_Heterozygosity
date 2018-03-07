#!/usr/bin/python
import sys
from fasta import readfasta

################################################################################
# Written by Homa Papoli - 14 Sept 2017
# Takes a single or a list of fasta sequences from a fasta file.
# Single fasta sequence:
# ./extract_seq_from_fasta.py file.fa single seqname > seqname.fa
# List of fasta sequences:
# ./extract_seq_from_fasta.py file.fa list seqlist > seqlist.fa
################################################################################

fasta = open(sys.argv[1],"r")
argument = sys.argv[2]

################################################################################
# Read the fasta file into a dictionary
################################################################################
fastaDict = readfasta(fasta)

# Function to output the fasta sequences into multiline file
def chunks(s, n):
	for start in range(0, len(s), n):
		yield s[start:start+n]

if argument == "single":
	sequence = sys.argv[3]
	print ">"+sequence
	for line in chunks(fastaDict[sequence], 50):
		print line	
else:
	sequence = open(sys.argv[3], "r")
	seqlist = []
	for line in sequence:
		line = line.strip("\n").replace(">","")
		seqlist.append(line)
	for element in seqlist:
		print ">"+element
		for line in chunks(fastaDict[element], 50):
			print line	
	
	
