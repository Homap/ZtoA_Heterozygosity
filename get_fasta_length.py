#!/usr/bin/python
import sys
from fasta import readfasta

#*******************************************************************************
# Written by Homa Papoli
# Reports the length of each fasta sequence.
# Usage: ./get_fasta_length file.fa > file.length.txt
#*******************************************************************************

# Open the fasta sequence
fasta = open(sys.argv[1], 'r')

# Read the fasta sequence into a dictionary
fastaseq = readfasta(fasta)

# Get the length of each sequence
for key in fastaseq.keys():
	print (key + '\t' + str(len(fastaseq.get(key))))
