#!/usr/bin/python
import sys
import operator
import numpy
from fasta import readfasta

#*********************************************************************************
# Written by Homa Papoli - 14 June 2016
# assembly_statistics.py takes a fasta file and their length. 
# It outputs length of all scaffolds, the number of scaffolds, name of the shortest
# and longest scaffolds, total length of the fasta sequences and N50.
# Usage: ./assembly_statistics.py file.fasta file.length.txt > file.statistics.txt
# file.fasta is the input fasta sequence.
# file.length.txt is the output fasta sequence length.
# file.statistics.txt is the output fasta statistics.
#*********************************************************************************

fasta = sys.argv[1]
length = open(sys.argv[2], 'w')

# Read the fasta sequence into a dictionary
with open(fasta, 'r') as f:
	fastaseq = readfasta(f)

# Get the length of each sequence and store it in length_dict
# key = scaffold_name, value = length of scaffold
length_dict = {}
for key in fastaseq.keys():
	if key in length_dict.keys():
		length_dict[key].append(len(fastaseq.get(key)))
	else:
		length_dict[key] = len(fastaseq.get(key))

# sort length_dict for the sequence length
sorted_length_dict_value = sorted(length_dict.items(), key = operator.itemgetter(1))

# store the scaffold length into length_list
length_list=[]
for element in sorted_length_dict_value:
	length_list.append(element[1])

# write scaffold and length in the output
for v,k in sorted_length_dict_value:
	length.write("%s\t%d\n" % (v, k))

# sort the list of value
length_list.sort()
# Print out some summary statistics
print ("Number of scaffolds: %d" % len(length_list))
print ("Shortest scaffold: %d" % min(length_list))
print ("Longest scaffold: %d is %s" % (max(length_list), max(length_dict.iteritems(), key=operator.itemgetter(1))[0]))
print ("Total length is: %d" % sum(length_list))
# Find N50 as described in http://www.broadinstitute.org/crd/wiki/index.php/N50
# go through the length list of scaffolds and add each unique element to the list unique
unique = []
for entry in length_list:
	if not entry in unique:
		unique.append(entry)
# for every number in the list unique, do as follows:
# [2, 2, 2, 3, 3, 4, 8, 8], multiply 3 by 2, 2 by 3, 1 by 4 and 2 by 8.
# then create a new list, n50 where there are 6 2s, 6 3s, 4 4s and 16 8s. 
n50 = []
for entry in unique:
	multiplier = length_list.count(entry) * entry
	for i in range(multiplier):
		n50.append(entry)

index = len(n50)/2
avg = []

# if index is an even number, then n50 would be the average of the index and one before it.
# if index is an odd number, n50 is the index then.
if index % 2==0:
	first = n50[index-1]
	second = n50[index]
	# avg.append(first)
	avg.append(second)
	n50 = numpy.mean(avg)
	print "The N50 is: %d" % n50
else:
	print "The N50 is: %d" % n50[index -1]

length.close()
