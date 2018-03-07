#!/usr/bin/python
from __future__ import division
import sys
from fasta import readfasta
from het import heterozygosity

fasta = sys.argv[1]

# Read fasta into a dictionary
with open(fasta, 'r') as f:
	fasta_dict = readfasta(f)

# Calculate heterozygosity for every scaffold
for key in fasta_dict.keys():
	het = heterozygosity(fasta_dict[key])
	print (key+"\t"+str(het[0])+"\t"+str(het[1])+"\t"+str(het[2]))




				

