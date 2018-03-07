#!/usr/bin/python
from __future__ import division
import sys
from het import heterozygosity
from fasta import readfasta
from intergenic import intergenicCoord
from SlidingWindow import slidingWindow

######################################################################
# The script uses 4 functions imported as modules:
# heterozyosity: function that return heterozygosity as the number of SNPs divided by the total sequenced length
# fasta: function that return the fasta sequence read as dictionary
# intergenic: function that return the intergenic sequence as a dictionary
# SlidingWindow: function that yield a list of sequence chunks cut from a the original fasta sequence based on window size and step size
######################################################################
# Usage: ./linked_selection.py Fasta.fa intergenic.txt winsize
# Formate of the files:
# Fasta.fa is a fasta file with sequences as single lines
# >Sequence_name
# ATGCATGCATGCWACGTNNN
# intergenic.txt is a text file with 4 columns. 1st column is the scaffold name, 2nd and 3rd columns are the start and end of the intergenic sequence and the 4th column is its length 
# scaffold_name start end length

fastaseq = open(sys.argv[1], "r") # File one is the fasta sequence
intergenic = open(sys.argv[2], "r") # File two is the intergenic file
winsize_i = int(sys.argv[3]) 

# Print header 
print ("Scaffold"+"\t"+"Orientation"+"\t"+"Window"+"\t"+"Window_Start"+"\t"+"Window_End"+"\t"+"Distance_from_Gene"+"\t"+"Intergenic_start"+"\t"+"Intergenic_end"+"\t"+"Num_SNPs"+"\t"+"Num_bases"+"\t"+"Het")

# Determine the window size and step size
winSize = winsize_i # How can I output this as variables?
step = winsize_i

# Read the fasta file into a dictionary
fastaDict = readfasta(fastaseq)

# Read the intergenic coordinates into a dictionary	
intergenDict = intergenicCoord(intergenic)

# Get a sequence in fasta that lies within the intergenic sequence -- V2
new_fastdict = {}
for scaffold in intergenDict.keys():
	counter  = 0
	if scaffold in fastaDict.keys():
		for coordinate in intergenDict.get(scaffold):
			counter = 0
			new_scaffold = fastaDict.get(scaffold)[int(coordinate[0]):int(coordinate[1])+1]
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
				print scaffold, "left", counter, start, end, distance, coordinate[0], coordinate[1], heterozygosity(i)
			counter = 0
			chunks_right = slidingWindow(sequence_right, winSize, step)
			for i in chunks_right:
				counter += 1
				end = str(int(coordinate[1]) - int(counter*winSize))
				start = str(int(end) + winSize)
				distance = str(int(coordinate[1]) - int(end))
				print scaffold, "right", counter, start, end, distance, coordinate[0], coordinate[1], heterozygosity(i)

fastaseq.close()
intergenic.close()
			

			



		

