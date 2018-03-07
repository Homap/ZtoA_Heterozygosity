#!/usr/bin/python
from __future__ import division
import sys
from het import heterozygosity
from fasta import readfasta
from SlidingWindow import slidingWindow

######################################################################
# The script uses 4 functions imported as modules:
# heterozyosity: function that return heterozygosity as the number of SNPs divided by the total sequenced length
# fasta: function that return the fasta sequence read as dictionary
# intergenic: function that return the intergenic sequence as a dictionary
# SlidingWindow: function that yield a list of sequence chunks cut from a the original fasta sequence based on window size and step size
######################################################################
# Usage: ./linked_selection.py Fasta.fa intergenic.txt
# Formate of the files:
# Fasta.fa is a fasta file with sequences as single lines
# >Sequence_name
# ATGCATGCATGCWACGTNNN
# intergenic.txt is a text file with 4 columns. 1st column is the scaffold name, 2nd and 3rd columns are the start and end of the intergenic sequence and the 4th column is its length 
# scaffold_name start end length

fastaseq = open(sys.argv[1], "r") # File one is the fasta sequence
intergenic = open(sys.argv[2], "r") # File two is the intergenic file

# Determine the window size and step size
winSize = int(sys.argv[3])
step = int(sys.argv[4])

# Read the fasta file into a dictionary
fastaDict = readfasta(fastaseq)
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

if winSize == 0 and step == 0:
	# Do not use sliding window approach, just calculate heterozygosity in the given region
	new_fastdict = {}
	for scaffold in intergenDict.keys():
		if scaffold in fastaDict.keys():
			for coordinate in intergenDict.get(scaffold):
				intron = fastaDict.get(scaffold)[int(coordinate[0]):int(coordinate[1])]
				print (scaffold+'\t'+str(coordinate[0])+'\t'+str(coordinate[1])+'\t'+str(heterozygosity(intron)[0])+'\t'+str(heterozygosity(intron)[1])+'\t'+str(heterozygosity(intron)[2]))
					

# Get a sequence in fasta that lies within the intergenic sequence -- V2
else:
	# Print header 
	print("Scaffold"+"\t"+"Orientation"+"\t"+"Window"+"\t"+"Window_Start"+"\t"+"Window_End"+"\t"+"Distance_from_Gene"+"\t"+"Intergenic_start"+"\t"+"Intergenic_end"+"\t"+"Num_SNPs"+"\t"+"Num_bases"+"\t"+"Het")
	new_fastdict = {}
	for scaffold in intergenDict.keys():
		counter  = 0
		if scaffold in fastaDict.keys():
#			print scaffold
			for coordinate in intergenDict.get(scaffold):
				counter = 0
				new_scaffold = fastaDict.get(scaffold)[int(coordinate[0]):int(coordinate[1])]
				length = len(new_scaffold)
				if length%2 == 1:
					length = length+1
				midpoint = int(length/2)
				sequence_left = new_scaffold[0:midpoint]
				sequence_right = new_scaffold[:(midpoint-1):-1]
				chunks_left = slidingWindow(sequence_left, winSize, step)
				for i in chunks_left:
					counter += 1
					end = str(int(coordinate[0]) + int(counter*winSize))
					start = str(int(end) - winSize)
					distance = str(abs(int(coordinate[0]) - int(end)))
					print (scaffold+'\t'+"left"+'\t'+str(counter)+'\t'+str(start)+'\t'+str(end)+'\t'+str(distance)+'\t'+str(coordinate[0])+'\t'+str(coordinate[1])+'\t'+str(heterozygosity(i)[0])+'\t'+str(heterozygosity(i)[1])+'\t'+str(heterozygosity(i)[2]))
				counter = 0
				chunks_right = slidingWindow(sequence_right, winSize, step)
				for i in chunks_right:
					counter += 1
					end = str(int(coordinate[1]) - int(counter*winSize))
					start = str(int(end) + winSize)
					distance = str(int(coordinate[1]) - int(end))
					print (scaffold+'\t'+"right"+'\t'+str(counter)+'\t'+str(start)+'\t'+str(end)+'\t'+str(distance)+'\t'+str(coordinate[0])+'\t'+str(coordinate[1])+'\t'+str(heterozygosity(i)[0])+'\t'+str(heterozygosity(i)[1])+'\t'+str(heterozygosity(i)[2]))

fastaseq.close()
intergenic.close()
			

			



		

