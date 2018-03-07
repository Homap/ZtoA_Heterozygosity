#!/usr/bin/python
import sys
from fasta import readfasta

# Written by Homa Papoli, Oct 2017
# Extract functional sequences from a fasta file and concatenates in order to 
# use them for resampling.
# Usage: ./compile_sequences.py seq.fasta seq.coordinates seq.scaf.chro 

fastaseq = open(sys.argv[1], "r") # File one is the fasta sequence
intergenic = open(sys.argv[2], "r") # File two is the intergenic file
chrfile = open(sys.argv[3], "r") # File three is the chromosome file that 
# contains all the scaffolds for that chromosome.

def chunks(s, n):
	for start in range(0, len(s), n):
		yield s[start:start+n]

# Read the fasta file into a dictionary
fastaDict = readfasta(fastaseq)
#print fastaDict
# Read the intergenic coordinates into a dictionary	
def intergenicCoord(intergenic):			
	intergenicDict= {}
	for line in intergenic:
		line = line.strip().split("\t")
		key, value = line[0], line[1:]
		if line[0] in intergenicDict.keys():
			intergenicDict[key].append(value)
		else:
			intergenicDict[key] = [value]
	return intergenicDict

intergenDict = intergenicCoord(intergenic)
#print intergenDict
# Read chromosome file into dictionary
chrdict = {}
for line in chrfile:
	line = line.strip("\n").split("\t")
	chromosome, scaffold = line[0], line[1]
	if chromosome in chrdict.keys():
		chrdict[chromosome].append(scaffold)
	else:
		chrdict[chromosome] = [scaffold]
#print intergenDict.keys()
# Extract the region of interest
new_fastdict = {}
for chromosome in chrdict.keys():
	l = list()
	for scaffold in chrdict[chromosome]:
#		print scaffold
	#for scaffold in intergenDict.keys():
		if scaffold in fastaDict.keys():
			if scaffold in intergenDict.keys(): # Not all scaffolds are annotated, so they're not found in the gff file.
#				print scaffold
				for coordinate in intergenDict[scaffold]:
					intron = fastaDict.get(scaffold)[int(coordinate[0]):int(coordinate[1])]
##					print intron
					l.append(intron)
##	print chromosome, l
	new_fastdict[chromosome] = ''.join(l)
#print new_fastdict

for chromosome in new_fastdict.keys():
	print ">"+chromosome
	for line in chunks(new_fastdict[chromosome],50):
		print line
		
		
fastaseq.close()
intergenic.close()
chrfile.close()

