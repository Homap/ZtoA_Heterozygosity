#!/usr/bin/python
import sys
from fasta import readfasta
from het import heterozygosity

# Usage:
# ./het_gene_density.py ../data/Cuculus_canorus.CDS_Repeats_Masked.fa ../data/cucu_out1_window > ../data/cucu_gene_density_het
# Fasta file must be masked for both repeats and CDS sequences. 
# 

fasta = sys.argv[1]
window = sys.argv[2]

# Read the fasta file into a dictionary
with open(fasta, 'r') as fastafile:
	fasta_dict = readfasta(fastafile)

# Read the window file into a dictionary
window_dict = {}
with open(window, 'r') as windowfile:
	for line in windowfile:
		line = line.strip('\n').split('\t')
		key, value = line[0], line[1:]
		if key in window_dict.keys():
			window_dict[key].append(value)
		else:
			window_dict[key] = [value]
# print window_dict

fasta_win = {}
for scaffold in window_dict.keys():
	for element in window_dict.get(scaffold):
		if scaffold in fasta_win.keys():
			fasta_win[scaffold].append(fasta_dict[scaffold][int(element[0])-1:int(element[1])])
		else:
			fasta_win[scaffold] = [fasta_dict[scaffold][int(element[0])-1:int(element[1])]]

# print fasta_win

for scaffold in fasta_win.keys():
	for sequence in fasta_win.get(scaffold):
		index = fasta_win[scaffold].index(sequence)
		print scaffold, window_dict[scaffold][index][0], window_dict[scaffold][index][1], window_dict[scaffold][index][2], heterozygosity(sequence)[0], heterozygosity(sequence)[1], heterozygosity(sequence)[2] 


