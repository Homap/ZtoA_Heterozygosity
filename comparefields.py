#!/usr/bin/python
# Created by: Homa Papoli
# Description: It compares two files based on two columns and sort one based on another.
#
# Usage: comparefields.py orderedfile unorderedfile output
#----------------------------------------------------------------------------------------
#========================================================================================
# imports
import sys
from datetime import datetime

# Store the current time:
startTime = datetime.now()

# Read the ordered file as a list
try:
	file1 = sys.argv[1]
	file2 = sys.argv[2]
	file3 = sys.argv[3]
	ordered = open(file1, 'r')
	unordered = open(file2, 'r')
	output = open(file3, 'w')
except IOError as ex:
	print (">> Sorry, Could not find the file: " + ex.stderror)

print ">> Reading the ordered file into a list"
orderedlist = ordered.readlines()
listoflist = []
for element in orderedlist:
	element = element.split('\t')
	combinedfields = element[0] + "_" + element[1]
	listoflist.append(combinedfields)
#print listoflist

print ">> Reading the unordered file into a dictionary"
unorderedlist = unordered.readlines()
unorddict = {}
for element in unorderedlist:
	element = element.split('\t')
	combinedfields = element[0] + "_" + element[1]
	unorddict[combinedfields] = element[2:]
#print unorddict

print ">> Writing the output"
for element in listoflist:
	if element in unorddict.keys():
		value = unorddict[element]
		element = element.split('_')
		output.write('\t'.join(element) + '\t' + '\t'.join(value).rstrip('\n') + '\n')
		
print ">> Done!"
ordered.close()
unordered.close()
output.close()

print datetime.now() - startTime
