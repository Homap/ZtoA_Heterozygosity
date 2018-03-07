#!/usr/bin/python
import sys
from fasta import readfasta

######################################################################################
# Written by Homa Papoli - 16 May 2016                                               #
# The script reads the CDS file produced by the ./extract_from_gff script into a     #
# dictionary, it then extracts the CDS coordinates from the fasta file and finally   #
# concatenates each sequence to the other to create a CDS with SNPs marked           #
# as IUPAC code.                                                                     #
# Usage: ./extract_CDS_from_Fasta.py cds.txt fasta.fa > output                       #  
######################################################################################

#*******************
# Specify the inputs
#*******************
cds_file = open(sys.argv[1], 'r')
fasta = open(sys.argv[2], 'r')

#************************************************************
# Read the fasta file into a dictionary
#************************************************************
fastaDict = readfasta(fasta)
# With the current readfasta() function, it's much faster to 
# use the single line fasta sequence
 
#*******************************************
# Read the CDS coordinates into a dictionary
#*******************************************
CDS_dict={}
for line in cds_file:
	line = line.strip('\n').split('\t')
	key, value = line[0], line[1:3]
	if key in CDS_dict.keys():
		CDS_dict[key].append(value)
	else:
		CDS_dict[key] = [value]
#print CDS_dict		
#*********************************************
# Filling the gene_dict with the CDS sequences
#*********************************************
gene_dict = {} # Define a dictionary called gene_dict
for keyseq in CDS_dict.keys(): # for each key in the CDS dictionary
	scaf = keyseq.split('_', 1)[0] # the key for each CDS coordinate appears as 'scaffold1_Cca_R000087', scaf = scaffold1
#	print scaf 
	if scaf in fastaDict.keys(): # if scaffold exists among the keys of fasta dictionary
		gene = '' # set gene to an empty string
		gene_name = keyseq.split('_', 1)[1] # set gene name to the second part of the key, in the example above, gene_name = Cca_R000087
		for element in CDS_dict[keyseq]: # for each value in the CDS dictionary
			if keyseq in gene_dict.keys(): # if the keyseq already exists in the gene_dict
				gene_dict[keyseq] = gene_dict[keyseq] + fastaDict[scaf][int(element[0]):int(element[1])] # add the sequence to the gene_name
			else: 
				gene_dict[keyseq] = fastaDict[scaf][int(element[0]):int(element[1])] 
			gene = gene + fastaDict[scaf][int(element[0]):int(element[1])]
#			print gene
						
#************************************************
# Print key, value of gene_dict in a fasta format
#************************************************
for key in gene_dict.keys():
	print '>' + key + '\n' + gene_dict[key]

cds_file.close()
fasta.close()
