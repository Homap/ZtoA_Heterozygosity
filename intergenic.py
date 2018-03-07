#!/usr/bin/python
import sys

# Function to read an intergenic coordinate file into a dictionary
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
