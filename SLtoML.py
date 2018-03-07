#!/usr/bin/python
import sys

################################################################################
# Written by Homa Papoli
# Converts a single line fasta file to multiline
# ./SLtoML.py infile.SL.fa outfile.ML.fa 
################################################################################
# Chunks function works as follows:
# Given a sequence: s = "NNNNNKAWQ", "for start in range(0, len(s), 2):", we
# have 0, 2, 4, 6, 8. This gives us: NN, NN, NK, AW, Q. 

def chunks(s, n):
	for start in range(0, len(s), n):
		yield s[start:start+n]

infile = open(sys.argv[1], 'r')
outfile = open(sys.argv[2], 'w')

for line in infile:
	line = line.rstrip('\n')
	if line.startswith('>'):
		outfile.write(line + '\n')
	else:
		for line in chunks(line, 50):
			outfile.write(line + '\n')		
