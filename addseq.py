#!/usr/bin/python

import sys

infile = open(sys.argv[1], 'r')
outfile = open(sys.argv[2], 'w')

seq_tail = 20000 * 'N'
for line in infile:
	if line.startswith('>'):
		outfile.write(line.rstrip('\n') + '\n')
	else:
		line = line.rstrip('\n').upper() + seq_tail
		outfile.write(line + '\n')
