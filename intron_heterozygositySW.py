#!/usr/bin/python
import sys
import collections
from fasta import readfasta
from het_SW import heterozygositySW

##########################################################################################
# Written by Homa Papoli - 31 May 2016                                               
# This script contains the following functions:
# readfasta()
# heterozygositySW()
# It reads the intron coordinates generated by extract_from_gff.py into a dictionary,
# extract the corresponding sequences from a fasta file and then runs the heterozygosity
# function to count the number of SNPs and bases.
# Usage: ./intron_heterozygositySW.py intron_coordinates.txt fasta.fa > output
##########################################################################################

#*************************************************** 
# Read the intron coordinates and the fasta sequence
#***************************************************
intron_coord = open(sys.argv[1], 'r')
fasta = open(sys.argv[2], 'r')

#******************************************
# Read the fasta sequence into a dictionary
#******************************************
fastaseq = readfasta(fasta)

#**********************************************
# Read the intron coordinates into a dictionary
#**********************************************
intron_dict = {}
for line in intron_coord:
	line = line.strip('\n').split('\t')
	key, value = line[0], line[1:]
	if key in intron_dict.keys():
		intron_dict[key].append(value)
	else:
		intron_dict[key] = [value]

#**************************************************************************		
# Extract the intronic sequences from fasta and read them into a dictionary
#**************************************************************************
for scaf in intron_dict.keys():
	scaffold = scaf.split('_', 1)[0]
	if scaffold in fastaseq.keys():
		for coordinate in intron_dict.get(scaf):
			# If the sequence is ABCD and the coordinates are 1 and 3, because the
			# indexing is from 0 in Python, we would take the sequence from 0 to 3
			# because the last base is exclusive.
			intron_seq = fastaseq.get(scaffold)[int(coordinate[0]):int(coordinate[1])]
			print (scaf+'\t'+str(coordinate[0])+'\t'+str(coordinate[1])+'\t'+str(heterozygositySW(intron_seq)[0])+'\t'+str(heterozygositySW(intron_seq)[1])+'\t'+str(heterozygositySW(intron_seq)[2]))

intron_coord.close()
fasta.close()
