#!/usr/bin/python
import sys
from fasta import readfasta

f = open(sys.argv[1], "r")

fd = readfasta(f)

for key in fd:
	sequence = fd[key][15720:15725]
	print sequence
	homozygote='ATCG'
	SNPs='RYSWKMBDHV'	
	homozygote_count = len([base.upper() for base in sequence if base.upper() in homozygote])		
	SNPs_count = len([base.upper() for base in sequence if base.upper() in SNPs])
	print homozygote_count
	print SNPs_count
