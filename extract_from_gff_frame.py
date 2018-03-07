#!/usr/bin/python
import sys
import collections

#########################################################################################################################
# Written by Homa Papoli - 26 May 2016                                                                                  #
# This script takes a gff file as an input and outputs the coordinates                                                  #
# for intergenic regions, intron or exons based on the term determined                                                  #
# in the argument.                               																	    #
# For intergenic, run:                                                  									    #
# Usage: ./extract_from_gff.py gff_file [intergenic/intron/CDS] [first/all/second_to_last] scaffold_length.txt > output #
# For CDS and intronic, run:																						    #
# Usage: ./extract_from_gff.py gff_file [intergenic/intron/CDS] [first/all/second_to_last] > output                     #
#########################################################################################################################

#*****************************************************************************
# Specify the inputs
#*****************************************************************************
gff = open(sys.argv[1], 'r')
element = sys.argv[2] # intergenic, intron or CDS
number = sys.argv[3] # if first intron/exon, enter 'first', if all introns/exons, enter 'all', 
					# if second to last introns/exons, enter 'second_to_last'
					# for intergenic, the third argument is always 'all'
					
# For the intergenic regions, the start and end of the scaffolds should also be
# taken into account. For this, we need the scaffold lengths.					
# Read scaffold lengths as a dictionary:
if element == "intergenic":
	scaf = open(sys.argv[4], "r")
	scaf_dict = {}
	for line in scaf:
		line = line.strip("\n").split("\t")
		scaffold, length = line[0], int(line[1])
		scaf_dict[scaffold] = length
	#print scaf_dict
#*****************************************************************************

#*****************************************************************************
# Read the mRNA and CDS coordinates into dictionaries
#*****************************************************************************
mRNA_coordinates = {}
CDS_coordinates = {}
for line in gff:
	line = line.strip('\n').split('\t')
	if line[2] == 'mRNA': # if the third column of the gff file is mRNA
		key, value = line[0], line[3:5] # set the scaffold name as key and the start:stop as value
		if key in mRNA_coordinates.keys(): # if the scaffold name is already present in the mRNA_coordinates dictionary
			mRNA_coordinates[key].append(value) # append the new coordinates to the list
		else:
			mRNA_coordinates[key] = [value] # otherwise set a new coordinate list to the new key
	elif line[2] == 'CDS': # if the third column of the gff file is CDS
		key, value = line[0]+'_'+line[8].replace(';', '').replace('Parent=', ''), [line[3], line[4], line[6], line[7]] # set the key in the form scaffold_genename and the start:stop as value 
		if key in CDS_coordinates.keys(): # if the third column of the gff file is CDS
			CDS_coordinates[key].append(value) 
		else:
			CDS_coordinates[key] = [value]
			
#print CDS_coordinates

#*****************************************************************************
# Read the intergenic coordinates into a dictionary
#*****************************************************************************			
#print mRNA_coordinates
if element == 'intergenic' and number == 'all': # element is the second argument of the command line, if it is intergenic
	intergenic = {} # define a dictionary called intergenic
	for scaffold in mRNA_coordinates.keys(): # look for the scaffold in the keys of mRNA_coordinates dictionary
		for i in range(len(mRNA_coordinates.get(scaffold))): # calculate the length of each value of the mRNA_coordinates dictionary and store its range in a variable called i 
			if i == len(mRNA_coordinates.get(scaffold)) - 1: # if i is equal to the length - 1 because i starts from 0 so if you have l = [1, 2, 3, 4], [i for i in range(len(l))] would be [0, 1, 2, 3]
				break # so the last element would have an index which is equal to the length of the list minus 1, if so, stop the loop
			else: # otherwise
				if scaffold in intergenic.keys(): # if scaffold exists among the keys in the intergenic dictionary
					intergenic[scaffold].append([mRNA_coordinates.get(scaffold)[i][1], mRNA_coordinates.get(scaffold)[i+1][0]]) # append to the list the intergenic coordinates
				else:
					intergenic[scaffold] = [[mRNA_coordinates.get(scaffold)[i][1], mRNA_coordinates.get(scaffold)[i+1][0]]]
					
	#**********************************************
	# Print the intergenic coordinate to the output
	#**********************************************
	# In order to obtain the intergenic sequences from the gff file, the sequences
	# between the two mRNA coordinates were obtained. However, the start and end
	# sequence of each scaffold is also potentially part of the intergenic sequence.
	# The code below takes into account these segments as well.			
	for key, value in intergenic.iteritems():
#		print key, value
#		Example when there are more than 2 coordinate sets:
#		scaffold3301 [['27574', '38089'], ['38253', '40600'], ['40764', '72871']]
#		scaffold3301    0       27573
#		scaffold3301    27573   38089
#		scaffold3301    38252   40600
#		scaffold3301    40763   72871
#		scaffold3301    72871   73321
		if len(value) > 2:
			for element in value:
				if value.index(element) == 0:
					print (key+'\t'+"0"+'\t'+str(int(element[0])-1))
					print (key+'\t'+str(int(element[0])-1)+'\t'+element[1])
				elif value.index(element) == len(value)-1:
					if key in scaf_dict.keys():
						print (key+'\t'+str(int(element[0])-1)+'\t'+element[1])
						print (key+'\t'+element[1]+'\t'+str(scaf_dict[key]))
				else:
					print (key+'\t'+str(int(element[0])-1)+'\t'+element[1])
#		scaffold54855 [['6377', '15720']]
#		scaffold54855   0       6376
#		scaffold54855   6376    15720
#		scaffold54855   15720   17673	
		elif len(value) == 1:
			for element in value:
				print (key+'\t'+"0"+'\t'+str(int(element[0])-1))
				print (key+'\t'+str(int(element[0])-1)+'\t'+element[1])
				if key in scaf_dict.keys():
					print (key+'\t'+element[1]+'\t'+str(scaf_dict[key]))
#		scaffold8907 [['56060', '68023'], ['85306', '102588']]
#		scaffold8907    0       56059
#		scaffold8907    56059   68023
#		scaffold8907    85305   102588
#		scaffold8907    102588  126423
		elif len(value) == 2:
			for element in value:
				if value.index(element) == 0:
					print (key+'\t'+"0"+'\t'+str(int(element[0])-1))
					print (key+'\t'+str(int(element[0])-1)+'\t'+element[1])
				elif value.index(element) == 1:
					if key in scaf_dict.keys():
						print (key+'\t'+str(int(element[0])-1)+'\t'+element[1])
						print (key+'\t'+element[1]+'\t'+str(scaf_dict[key]))
			
#*****************************************************************************
# Read the intronic coordinates into a dictionary
#*****************************************************************************		
if element == 'intron':	# if element is intron
	intronic = {} # define a dictionary called intron
	for scaffold in CDS_coordinates.keys(): # look for the scaffold in the keys of CDS_coordinates dictionary
		if len(CDS_coordinates[scaffold]) != 1: # if length of the value list is not equal to 1, because if it is equal to 1, it means there is only one CDS and no intron
			for position1, position2 in zip(CDS_coordinates[scaffold], CDS_coordinates[scaffold][1:]): # now use the zip function which pairs every two elements in a list
				if scaffold in intronic.keys():
					intronic[scaffold].append([position1[1], position2[0]])
				else:	
					intronic[scaffold] = [[position1[1], position2[0]]]
					
	#*********************************************
	# Print all intronic coordinates to the output
	#*********************************************				
	if number == 'all': # in case the second argument is intron, the third argument can take 'all', 'first' or 'second_to_last'
		for key, value in intronic.iteritems():
			for element in value:
				print (key+'\t'+str(int(element[0])-1)+'\t'+element[1])
				
	#***************************************************
	# Print the first intronic coordinates to the output
	#***************************************************
	if number == 'first': # if it is first, it only ouputs the coordinates of the first intron
		for key, value in intronic.iteritems():
			if len(value)>1:
				print (key+'\t'+str(int(value[0][0])-1)+'\t'+value[0][1])
				
	#************************************************************
	# Print from the second to the last coordinates to the output
	#************************************************************				
	if number == 'second_to_last': # if it is second_to_last it outputs the intronic coordinates from the second to the last one
		for key, value in intronic.iteritems():
			if len(value)>1:
				for index, item in enumerate(value):
					if index != 0:
						print (key+'\t'+str(int(item[0])-1)+'\t'+item[1])
			elif len(value) == 1:
				print (key+'\t'+str(int(value[0][0])-1)+'\t'+value[0][1])
						
#************************************************
# Print all the CDS coordinates into a dictionary
#************************************************							
if element == 'CDS' and number =='all': # the same thing applies to the CDS, if it is all, it prints the coordinates for all CDS
	for key, value in CDS_coordinates.iteritems():
		if value[0][2] == "+": # if CDS is in forward strand
			for index, item in enumerate(value): # loop over the list
				if index == 0: # if index is equal to 0, that is we are at the first element of the list
					if value[index][3] == "0": # if first base of CDS is the first base of the codon
						print (key+'\t'+str(int(value[index][0])-1)+'\t'+value[index][1]+'\t'+value[index][2]+'\t'+value[index][3]) # gff is 1 base, so, if in gff, CDS start is at 24190 and end at 24271, to extract the sequence in python, I substract 1 from start but keep the end as it is because in python, we have [24190, 24271), that is closed start and open end of a set.
					elif value[index][3] == "1": # if second base of CDS is the first base of the codon
						print (key+'\t'+str(value[index][0])+'\t'+value[index][1]+'\t'+value[index][2]+'\t'+value[index][3])
					elif value[index][3] == "2": # if third base of CDS is the first base of the codon
						print (key+'\t'+str(int(value[index][0])+1)+'\t'+value[index][1]+'\t'+value[index][2]+'\t'+value[index][3])
				else: # if index is not equal to 0
					print (key+'\t'+str(int(value[index][0])-1)+'\t'+value[index][1]+'\t'+value[index][2]+'\t'+value[index][3])
		elif value[0][2] == "-": # if CDS is in reverse strand
			for index, item in enumerate(value):
				if index == len(value)-1:
					if value[index][3] == "0": 
						print (key+'\t'+str(int(value[index][0])-1)+'\t'+value[index][1]+'\t'+value[index][2]+'\t'+value[index][3])
					elif value[index][3] == "1": # if first base of CDS is the first base of the codon
						print (key+'\t'+str(int(value[index][0])-1)+'\t'+str(int(value[index][1])-1)+'\t'+value[index][2]+'\t'+value[index][3])
					elif value[index][3] == "2":
						print (key+'\t'+str(int(value[index][0])-1)+'\t'+str(int(value[index][1])-2)+'\t'+value[index][2]+'\t'+value[index][3])
				else: # if index is not equal to the index of the last element
					print (key+'\t'+str(int(value[index][0])-1)+'\t'+value[index][1]+'\t'+value[index][2]+'\t'+value[index][3])
		
#*************************************************
# Print the first CDS coordinate into a dictionary
#*************************************************		
elif element == 'CDS' and number == 'first': # if it is first, it prints the coordinates for the first CDS only
	for key, value in CDS_coordinates.iteritems():
		if len(value)>1:
			print (key+'\t'+str(int(value[0][0])-1)+'\t'+value[0][1])
			
#**************************************************************
# Print the second to the last CDS coordinate into a dictionary
#**************************************************************
elif element == 'CDS' and number == 'second_to_last': # if it is second_to_last, it prints the CDS coordinates for the second one to the last
	for key, value in CDS_coordinates.iteritems():
		if len(value)>1:
			for index, item in enumerate(value):
				if index !=0:
					print (key+'\t'+str(int(item[0])-1)+'\t'+item[1])
		elif len(value) == 1:
			print (key+'\t'+str(int(value[0][0])-1)+'\t'+value[0][1])
	
			
gff.close()
if element == "intergenic":
	scaf.close()
