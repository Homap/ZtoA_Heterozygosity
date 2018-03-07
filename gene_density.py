#!/usr/bin/python
import sys

# Usage: ./gene_density.py species_scaffold_length winsize > species_window

scaflength= open(sys.argv[1], 'r')

winSize = int(sys.argv[2])

for line in scaflength:
	line = line.strip('\n').split('\t')
	scaflen = int(line[1])
	# if length of the scaffold is smaller than the window size, ignore that scaffold.
	if scaflen < winSize:
		pass
	else:
	 	if scaflen%winSize == 0:
			for chunk in range(0, scaflen, winSize):
				print (str(line[0])+"\t"+str(chunk+1)+"\t"+str(chunk+winSize))
		elif scaflen%winSize != 0:
			for chunk in range(0, scaflen, winSize):
				if chunk+winSize <= scaflen:
					print (str(line[0])+"\t"+str(chunk+1)+"\t"+str(chunk+winSize))
				else:
					print (str(line[0])+"\t"+str(chunk+1)+"\t"+str(chunk+scaflen%winSize))
				
scaflength.close()
