#!/usr/bin/python
import sys
#import random
import numpy as np
from fasta import readfasta
from het import heterozygosity

#*******************************************************************************
# Written by Homa Papoli - October 2017
#*******************************************************************************
# Script contains functions to:
# 1. Generate a fasta sequence without N
# 2. Perform resampling from the fasta sequence to generate new sequence
# 3. Calculate heterozygosity from the new sequence
#*******************************************************************************

f1 = open(sys.argv[1], "r")
seq = sys.argv[2] # indicate which chromosome to resample
replicates = int(sys.argv[3]) # number of resampling
num = int(sys.argv[4]) # indicate the length from which to sample

# Read the fasta file into a dictionary. 
fastadict = readfasta(f1)
#def resampling_f(fastadict, seq, num):
#	fastadict[seq] = fastadict[seq].replace("N","").replace("n","")
#	l = []
#	# If sampling the sequence as long as the original one	
#	# new_seq = ''.join([random.choice(fastadict[seq]) for nuc in fastadict[seq]])
#	# If sampling the sequence for a specific set of number
#	new_seq = ''.join([random.choice(fastadict[seq]) for i in range(num)]) # New sequences
#	new_seq_het = list(heterozygosity(new_seq))[2] # Het of the new sequence
#	l.append(new_seq_het)
#	return l

def resampling_f(fastadict, seq, n, k):
	fastadict[seq] = fastadict[seq].replace("N","").replace("n","")
	seq_list = np.random.choice(tuple(fastadict[seq]), replace=True, size=(n * k,)).view('S{k}'.format(k=k))
	l = []
	for element in seq_list:
		new_seq_het = list(heterozygosity(element))[2] # Het of the new sequence
		l.append(new_seq_het)
	return l

for i in resampling_f(fastadict, seq, replicates, num):
	print i
	

	
