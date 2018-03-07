#!/usr/bin/python
import sys

###########################################################
# Written by Homa Papoli in March 2016
# This script gets a mRNA coordinates file and outputs the coordinates for the intergenic regions.
# Usage: ./mRNA_to_intergenic.py mRNA.txt intergenic.txt
# mRNA.txt is produced from gff files by greping for 'mRNA'

mRNA = open(sys.argv[1], 'r')
intergenic = open(sys.argv[2], 'w')

CDS_coordinates={}
for line in mRNA:
	line = line.strip().split("\t")
	key, value = line[0], line[3:5]
	if line[0] in CDS_coordinates.keys():
		CDS_coordinates[key].append(value)
	else:
		CDS_coordinates[key] = [value]

for scaffold in CDS_coordinates.keys():
	for i in range(len(CDS_coordinates.get(scaffold))):
		if i == len(CDS_coordinates.get(scaffold)) - 1:
			break
		else:
			intergenic.write(scaffold + "\t" + str(CDS_coordinates.get(scaffold)[i][1]) + "\t" + str(CDS_coordinates.get(scaffold)[i+1][0]) + "\t" + str(abs(int(CDS_coordinates.get(scaffold)[i+1][0]) - int(CDS_coordinates.get(scaffold)[i][1]))) +  "\n")

mRNA.close()	
