#!/usr/bin/python
from __future__ import division
import sys
import collections
import warnings
from fasta import readfasta

#*******************
# Specify the inputs
#*******************
#inFile = open(sys.argv[1], 'r')

#************************************************************
# Read the fasta file into a dictionary
#************************************************************
#fastaseq = readfasta(inFile)
#************************************************************    
# Specify four nucleotides
#************************************************************
nucs = ["A", "T", "C", "G"]
#************************************************************
# Specify IUPAC codes as a dictionary
#************************************************************
IUPAC_code = {'R':['A', 'G'], 'Y':['C', 'T'], 'S':['G', 'C'], 'W':['A', 'T'], 'K':['G', 'T'], 'M':['A', 'C']}#, 'B':['C', 'G', 'T'], 'D':['A', 'G', 'T'], 'H':['A', 'C', 'T'], 'V':['A', 'C', 'G']}
#************************************************************
# Four fold degenerate sites
four_fold_codons = {
'A': ['GCT', 'GCC', 'GCA', 'GCG'],
'G': ['GGT', 'GGC', 'GGA', 'GGG'],  
'P': ['CCT', 'CCC', 'CCA', 'CCG'], 
'T': ['ACT', 'ACC', 'ACA', 'ACG'], 
'V': ['GTT', 'GTC', 'GTA', 'GTG'],
'R': ['CGT', 'CGC', 'CGA', 'CGG'],
'L': ['CTT', 'CTC', 'CTA', 'CTG'],
'S': ['TCT', 'TCC', 'TCA', 'TCG']
}
#************************************************************
# Specify amino acid codons in a dictionary
#************************************************************
DNA_CODON_TABLE = {
	'AAA':'K', 'AAC':'N', 'AAG':'K', 'AAR':'K', 'AAT':'N', 'AAY':'N', 'ACA':'T', 'ACB':'T', 
	'ACC':'T', 'ACD':'T', 'ACG':'T', 'ACH':'T', 'ACK':'T', 'ACM':'T', 'ACN':'T', 'ACR':'T', 
	'ACS':'T', 'ACT':'T', 'ACV':'T', 'ACW':'T', 'ACY':'T', 'AGA':'R', 'AGC':'S', 'AGG':'R', 
	'AGR':'R', 'AGT':'S', 'AGY':'S', 'ATA':'I', 'ATC':'I', 'ATG':'M', 'ATH':'I', 'ATM':'I', 
	'ATT':'I', 'ATW':'I', 'ATY':'I', 'CAA':'Q', 'CAC':'H', 'CAG':'Q', 'CAR':'Q', 'CAT':'H', 
	'CAY':'H', 'CCA':'P', 'CCB':'P', 'CCC':'P', 'CCD':'P', 'CCG':'P', 'CCH':'P', 'CCK':'P', 
	'CCM':'P', 'CCN':'P', 'CCR':'P', 'CCS':'P', 'CCT':'P', 'CCV':'P', 'CCW':'P', 'CCY':'P', 
	'CGA':'R', 'CGB':'R', 'CGC':'R', 'CGD':'R', 'CGG':'R', 'CGH':'R', 'CGK':'R', 'CGM':'R', 
	'CGN':'R', 'CGR':'R', 'CGS':'R', 'CGT':'R', 'CGV':'R', 'CGW':'R', 'CGY':'R', 'CTA':'L', 
	'CTB':'L', 'CTC':'L', 'CTD':'L', 'CTG':'L', 'CTH':'L', 'CTK':'L', 'CTM':'L', 'CTN':'L', 
	'CTR':'L', 'CTS':'L', 'CTT':'L', 'CTV':'L', 'CTW':'L', 'CTY':'L', 'GAA':'E', 'GAC':'D', 
	'GAG':'E', 'GAR':'E', 'GAT':'D', 'GAY':'D', 'GCA':'A', 'GCB':'A', 'GCC':'A', 'GCD':'A', 
	'GCG':'A', 'GCH':'A', 'GCK':'A', 'GCM':'A', 'GCN':'A', 'GCR':'A', 'GCS':'A', 'GCT':'A', 
	'GCV':'A', 'GCW':'A', 'GCY':'A', 'GGA':'G', 'GGB':'G', 'GGC':'G', 'GGD':'G', 'GGG':'G', 
	'GGH':'G', 'GGK':'G', 'GGM':'G', 'GGN':'G', 'GGR':'G', 'GGS':'G', 'GGT':'G', 'GGV':'G', 
	'GGW':'G', 'GGY':'G', 'GTA':'V', 'GTB':'V', 'GTC':'V', 'GTD':'V', 'GTG':'V', 'GTH':'V', 
	'GTK':'V', 'GTM':'V', 'GTN':'V', 'GTR':'V', 'GTS':'V', 'GTT':'V', 'GTV':'V', 'GTW':'V', 
	'GTY':'V', 'MGA':'R', 'MGG':'R', 'MGR':'R', 'TAA':'_', 'TAC':'Y', 'TAG':'_', 'TAR':'_', 
	'TAT':'Y', 'TAY':'Y', 'TCA':'S', 'TCB':'S', 'TCC':'S', 'TCD':'S', 'TCG':'S', 'TCH':'S', 
	'TCK':'S', 'TCM':'S', 'TCN':'S', 'TCR':'S', 'TCS':'S', 'TCT':'S', 'TCV':'S', 'TCW':'S', 
	'TCY':'S', 'TGA':'_', 'TGC':'C', 'TGG':'W', 'TGT':'C', 'TGY':'C', 'TRA':'_', 'TTA':'L', 
	'TTC':'F', 'TTG':'L', 'TTR':'L', 'TTT':'F', 'TTY':'F', 'YTA':'L', 'YTG':'L', 'YTR':'L',
	'---': '-', '...': '-', '~~~': '-'
}
#*******************************************************************************
# Stop codons
#*******************************************************************************
stop_codons = ['TGA', 'TAA', 'TAG', 'TRA', 'TAR']
#*******************************************************************************



# Get codons from a DNA sequence
#*******************************************************************************
def get_codons(sequence, nucs=nucs, DNA_CODON_TABLE=DNA_CODON_TABLE):
	all_codons = []
	for i in range(0, len(sequence)-2, 3): # a codon every three bases
		codon = sequence[i:i+3]
		if not 'N' in codon:
			all_codons.append(codon) # append all codons in the all_codons list
	return all_codons
#*******************************************************************************
# Count synonymous-nonsynonymous sites
#*******************************************************************************
def count_sites(sequence, nucs=nucs, DNA_CODON_TABLE=DNA_CODON_TABLE):
	codons = get_codons(sequence) # get_codons returns a list of all codons
	syn_sitesl = [] # list of synonymous sites
	nonsyn_sitesl = []  # list of non-synonymous sites
	for codon in codons: # for each codon in codons list
		# Synonymous sites with SNP ********************************************
		matching = [i for i in codon if i in IUPAC_code.keys()] # check if codon contains a SNP
		if matching: # if there is a SNP in codon
			if len(matching) == 1: # if there is only one SNP in a codon
				# example: AGY would be AGC and AGT
				syn_nonsyn1 = [] # append "s" or "ns" for the first codon = AGC
				syn_nonsyn2 = [] # append "s" or "ns" for the second codon = AGT
				# **************************************************************
				SNP = list(matching)[0] # the IUPAC code for the SNP in the codon
				codon1 = codon.replace(SNP, IUPAC_code[SNP][0]) # replace the IUPAC code, e.g. Y with C
				aa1 = DNA_CODON_TABLE[codon1] # get the respective amino acid from codon table
				codon2 = codon.replace(SNP, IUPAC_code[SNP][1]) # replace the IUPAC code, e.g. Y with T
				aa2 = DNA_CODON_TABLE[codon2] # get the respective amino acid from codon table
				for index in range(0,3,1): # each codon has 3 bases, loop over them 
					for nuc in nucs: # loop over the 4 nucleotides
						if not codon1[index] == nuc: # if the base in codon is not the base in the nucs
							codon1_l = list(codon1) # set codon1 as a list
							codon1_l[index] = nuc # change nucleotide of the codon in a given position
							new_codon1 = ''.join(codon1_l) # change the codon list back into string
							if not new_codon1 in stop_codons: # if the new codon is not in the stop codons
								if DNA_CODON_TABLE[new_codon1] == aa1: # check if the changed codon would translate into the same amino acid
									syn_nonsyn1.append("s") # in that case append "s" it to syn_nonsyn1 list
								else: # if the amino acid changes
									syn_nonsyn1.append("ns") # append "ns" to syn_nonsyn1 list 
						if not codon2[index] == nuc:
							codon2_l = list(codon2)
							codon2_l[index] = nuc
							new_codon2 = ''.join(codon2_l)
							if not new_codon2 in stop_codons:
								if DNA_CODON_TABLE[new_codon2] == aa2:
									syn_nonsyn2.append("s")
								else:
									syn_nonsyn2.append("ns")
				# do calculations for synonymous and nonsynonymous sites
				syn_sites1 = 3*(syn_nonsyn1.count("s")/len(syn_nonsyn1))
				nonsyn_sites1 = 3*(syn_nonsyn1.count("ns")/len(syn_nonsyn1))
				syn_sites2 = 3*(syn_nonsyn2.count("s")/len(syn_nonsyn2))
				nonsyn_sites2 = 3*(syn_nonsyn2.count("ns")/len(syn_nonsyn2))
				SNP_syn_sites = (syn_sites1+syn_sites2)/2
				SNP_nonsyn_sites = (nonsyn_sites1+nonsyn_sites2)/2
				syn_sitesl.append(SNP_syn_sites)
				nonsyn_sitesl.append(SNP_nonsyn_sites)
		else: # if the codon does not contain a SNP
			syn_nonsyn = []
			aa = DNA_CODON_TABLE[codon] # get the amino acid of the codon
			for index in range(0,3,1): # for every base in the codon					
				for nuc in nucs: # for every base in A, T, G, C
					if not codon[index] == nuc:
						codon_l = list(codon)
						codon_l[index] = nuc
						new_codon = ''.join(codon_l)
						if not new_codon in stop_codons:
							if DNA_CODON_TABLE[new_codon] == aa:
								syn_nonsyn.append("s")
							else:
								syn_nonsyn.append("ns")
			syn_sites = 3*(syn_nonsyn.count("s")/len(syn_nonsyn))
			nonsyn_sites = 3*(syn_nonsyn.count("ns")/len(syn_nonsyn))
			syn_sitesl.append(syn_sites)
			nonsyn_sitesl.append(nonsyn_sites)
	return sum(syn_sitesl), sum(nonsyn_sitesl)
#*******************************************************************************	
# Count four-fold sites
#*******************************************************************************		
def count_fourfold_sites (sequence, nucs=nucs, DNA_CODON_TABLE=DNA_CODON_TABLE):
	codons = get_codons(sequence) # get_codons returns a list of all codons
	four_fold = [] #list of four-fold degenerate sites
	for codon in codons: # for each codon in codons list
		matching = [i for i in codon if i in IUPAC_code.keys()]  #check if codon contains a SNP
		if matching:  #if there is a SNP in codon
			if len(matching) == 1:  #if there is only one SNP in a codon
				SNP = list(matching)[0]
				codon1 = codon.replace(SNP, IUPAC_code[SNP][0])
				aa1 = DNA_CODON_TABLE[codon1]
				codon2 = codon.replace(SNP, IUPAC_code[SNP][1])
				aa2 = DNA_CODON_TABLE[codon2]
				# below, it is not enough to check if an amino acid is among the keys of the four_fold_codons.
				# two amino acids, Arg(L) and Leu(L) have 6 codons. One must check if a given codon is among the 4 listed in four_fold_codons
				# as the value of L and R. 
				if aa1 == aa2: # if two amino acids from two codons are the same
					if aa1 in four_fold_codons.keys():
						if codon1 in four_fold_codons[aa1]: # if one of the codons is among the codons of the amino acid
							four_fold.append(1) # there exists a four fold degenerate site
				else:
					if aa1 in four_fold_codons.keys() and aa2 in four_fold_codons.keys(): # if two amino acids are different but each is an amino acid with four codons
						if codon1 in four_fold_codons[aa1] and codon2 in four_fold_codons[aa2]: # if both codons are found in the codons of the amino acid
							four_fold.append(1) # there exists a four fold degenerate site
						elif codon1 in four_fold_codons[aa1] or codon2 in four_fold_codons[aa2]: # if only one of the codons is found 
							four_fold.append(0.5) # there exists 0.5 four fold degenerate site
						else:
							four_fold.append(0) # otherwise, there is no four fold degenerate site
					elif aa1 in four_fold_codons.keys() or aa2 in four_fold_codons.keys(): # if only one of the amino acids is in four_fold_codons
						if codon1 in four_fold_codons.get(aa1,"X") or codon2 in four_fold_codons.get(aa2,"X"): # if at least one of them contains a codon that is in four_fold_codons
							four_fold.append(0.5) # there exists 0.5 four fold degenerate site
		else: # if there is no SNP in the codon
			aa = DNA_CODON_TABLE[codon]	# take the amino acid 									
			if aa in four_fold_codons.keys(): # if amino acid is among the keys of the four_fold_codons
				if codon in four_fold_codons[aa]: # if codon is among the codons of the amino acid found in four_fold_codons
					four_fold.append(1) # there exists 1 four fold degenerate site
	return sum(four_fold)
#*******************************************************************************
# Count zero-fold sites
#*******************************************************************************	
def count_zerofold_sites (sequence, nucs=nucs, DNA_CODON_TABLE=DNA_CODON_TABLE):
	codons = get_codons(sequence) # get_codons returns a list of all codons
	zero_fold = [] #list of zero-fold degenerate sites
	for codon in codons:
		if codon in stop_codons:
			if not codon in codons[-1]: #raise error if there is a stop codon in the middle of the ORF
				warnings.warn("Truncated ORF: Stop codon found in the middle of ORF!")
			elif codon in codons[-1]: # pass if there is a stop codon at the end of the ORF
				pass
		else:
			matching = [i for i in codon if i in IUPAC_code.keys()] #check if codon contains a SNP
			if matching:  #if there is a SNP in codon
				if len(matching) == 1: #if there is only one SNP in a codon
					SNP = list(matching)[0]
					codon1 = codon.replace(SNP, IUPAC_code[SNP][0])
					aa1 = DNA_CODON_TABLE[codon1]
					codon2 = codon.replace(SNP, IUPAC_code[SNP][1])
					aa2 = DNA_CODON_TABLE[codon2]
#					0-fold degenerate sites with SNP
					for index in range(0,3,1): # loop over the 3 bases in the codon
						zero_fold1 = [] # first list for zero fold 
						zero_fold2 = [] # second list for zero fold
						for nuc in nucs: # loop over A, C, G, T
							if not codon1[index] == nuc:
								codon1_l = list(codon1)
								codon1_l[index] = nuc
								new_codon1 = ''.join(codon1_l)
								if not new_codon1 in stop_codons: # if the new_codon1 is not in stop codons
									if DNA_CODON_TABLE[new_codon1] == aa1:
										zero_fold1.append("s") # add "s" in zero_fold1 
									else: 
										zero_fold1.append("ns") # add "ns" in zero_fold2
							if not codon2[index] == nuc: 
								codon2_l = list(codon2)
								codon2_l[index] = nuc
								new_codon2 = ''.join(codon2_l)
								if not new_codon2 in stop_codons: # if the new_codon2 is not in stop codons
									if DNA_CODON_TABLE[new_codon2] == aa2:
										zero_fold2.append("s")
									else:
										zero_fold2.append("ns")
						if list(set(zero_fold1))[0] == "ns" and list(set(zero_fold2))[0] == "ns":
							zero_fold.append(1)
						elif list(set(zero_fold1))[0] == "ns" or list(set(zero_fold2))[0] == "ns":
							zero_fold.append(0.5)
#					0-fold degenerate sites no SNP
			else:			
				aa = DNA_CODON_TABLE[codon]
				for index in range(0,3,1):
					zero_fold_per_codon = []
					for nuc in nucs:
						if not codon[index] == nuc:
							codon_l = list(codon)
							codon_l[index] = nuc
							new_codon = ''.join(codon_l)
							if not new_codon in stop_codons:
								if DNA_CODON_TABLE[new_codon] == aa:
									zero_fold_per_codon.append("s")
								else:
									zero_fold_per_codon.append("ns")
					if list(set(zero_fold_per_codon))[0] == "ns":
						zero_fold.append(1)
	return sum(zero_fold)
#*******************************************************************************
# Count synonymous-nonsynonymous sites - SW
#*******************************************************************************
def count_sites_SW(sequence, nucs=nucs, DNA_CODON_TABLE=DNA_CODON_TABLE):
	codons = get_codons(sequence) # get_codons returns a list of all codons
	syn_sitesl = [] # list of synonymous sites
	nonsyn_sitesl = []  # list of non-synonymous sites
	for codon in codons: # for each codon in codons list
#		print codon
		# Synonymous sites with SNP ********************************************
		matching = [i for i in codon if i in IUPAC_code.keys()] # check if codon contains a SNP
		if matching: # if there is a SNP in codon
			if len(matching) == 1: # if there is only one SNP in a codon
				# example: AGY would be AGC and AGT
				syn_nonsyn1 = [] # append "s" or "ns" for the first codon = AGC
				syn_nonsyn2 = [] # append "s" or "ns" for the second codon = AGT
				# **************************************************************
				SNP = list(matching)[0] # the IUPAC code for the SNP in the codon
				codon1 = codon.replace(SNP, IUPAC_code[SNP][0]) # replace the IUPAC code, e.g. Y with C
#				print codon1
				aa1 = DNA_CODON_TABLE[codon1] # get the respective amino acid from codon table
#				print aa1
				codon2 = codon.replace(SNP, IUPAC_code[SNP][1]) # replace the IUPAC code, e.g. Y with T
#				print codon2
				aa2 = DNA_CODON_TABLE[codon2] # get the respective amino acid from codon table
#				print aa2
				for index in range(0,3,1): # each codon has 3 bases, loop over them 
				# Codon 1
					if codon1[index] == "A":
						codon1_l = list(codon1)
						codon1_l[index] = "T"
						new_codon1 = ''.join(codon1_l)
#						print new_codon1
						if not new_codon1 in stop_codons:
							if DNA_CODON_TABLE[new_codon1] == aa1:
								syn_nonsyn1.append("s")
							else:
								syn_nonsyn1.append("ns") 
#								print DNA_CODON_TABLE[new_codon1]
					elif codon1[index] == "T":
						codon1_l = list(codon1)
						codon1_l[index] = "A"
						new_codon1 = ''.join(codon1_l)
						if not new_codon1 in stop_codons:
							if DNA_CODON_TABLE[new_codon1] == aa1:
								syn_nonsyn1.append("s")
							else:
								syn_nonsyn1.append("ns") 								
					elif codon1[index] == "C":
						codon1_l = list(codon1)
						codon1_l[index] = "G"
						new_codon1 = ''.join(codon1_l)
						if not new_codon1 in stop_codons:
							if DNA_CODON_TABLE[new_codon1] == aa1:
								syn_nonsyn1.append("s")
							else:
								syn_nonsyn1.append("ns") 								
					elif codon1[index] == "G":
						codon1_l = list(codon1)
						codon1_l[index] = "C"
						new_codon1 = ''.join(codon1_l)
						if not new_codon1 in stop_codons:
							if DNA_CODON_TABLE[new_codon1] == aa1:
								syn_nonsyn1.append("s")
							else:
								syn_nonsyn1.append("ns") 								
					# Codon 2
					if codon2[index] == "A":
						codon2_l = list(codon2)
						codon2_l[index] = "T"
						new_codon2 = ''.join(codon2_l)
						if not new_codon2 in stop_codons:
							if DNA_CODON_TABLE[new_codon2] == aa2:
								syn_nonsyn2.append("s")
							else:
								syn_nonsyn2.append("ns") 
					elif codon2[index] == "T":
						codon2_l = list(codon2)
						codon2_l[index] = "A"
						new_codon2 = ''.join(codon2_l)
						if not new_codon2 in stop_codons:
							if DNA_CODON_TABLE[new_codon2] == aa2:
								syn_nonsyn2.append("s")
							else:
								syn_nonsyn2.append("ns") 								
					elif codon2[index] == "C":
						codon2_l = list(codon2)
						codon2_l[index] = "G"
						new_codon2 = ''.join(codon2_l)
						if not new_codon2 in stop_codons:
							if DNA_CODON_TABLE[new_codon2] == aa2:
								syn_nonsyn2.append("s")
							else:
								syn_nonsyn2.append("ns") 								
					elif codon2[index] == "G":
						codon2_l = list(codon2)
						codon2_l[index] = "C"
						new_codon2 = ''.join(codon2_l)
						if not new_codon2 in stop_codons:
							if DNA_CODON_TABLE[new_codon2] == aa2:
								syn_nonsyn2.append("s")
							else:
								syn_nonsyn2.append("ns")
				# do calculations for synonymous and nonsynonymous sites
				syn_sites1 = 3*(syn_nonsyn1.count("s")/len(syn_nonsyn1))
				nonsyn_sites1 = 3*(syn_nonsyn1.count("ns")/len(syn_nonsyn1))
				syn_sites2 = 3*(syn_nonsyn2.count("s")/len(syn_nonsyn2))
				nonsyn_sites2 = 3*(syn_nonsyn2.count("ns")/len(syn_nonsyn2))
				SNP_syn_sites = (syn_sites1+syn_sites2)/2
				SNP_nonsyn_sites = (nonsyn_sites1+nonsyn_sites2)/2
				syn_sitesl.append(SNP_syn_sites)
				nonsyn_sitesl.append(SNP_nonsyn_sites)
		else: # if the codon does not contain a SNP
			syn_nonsyn = []
			aa = DNA_CODON_TABLE[codon] # get the amino acid of the codon
			for index in range(0,3,1): # for every base in the codon
				if codon[index] == "A":
					codon_l = list(codon)
					codon_l[index] = "T"
					new_codon = ''.join(codon_l)
					if not new_codon in stop_codons:
						if DNA_CODON_TABLE[new_codon] == aa:
							syn_nonsyn.append("s")
						else:
							syn_nonsyn.append("ns")
				if codon[index] == "T":
					codon_l = list(codon)
					codon_l[index] = "A"
					new_codon = ''.join(codon_l)
					if not new_codon in stop_codons:
						if DNA_CODON_TABLE[new_codon] == aa:
							syn_nonsyn.append("s")
						else:
							syn_nonsyn.append("ns")
				if codon[index] == "C":
					codon_l = list(codon)
					codon_l[index] = "G"
					new_codon = ''.join(codon_l)
					if not new_codon in stop_codons:
						if DNA_CODON_TABLE[new_codon] == aa:
							syn_nonsyn.append("s")
						else:
							syn_nonsyn.append("ns")
				if codon[index] == "G":
					codon_l = list(codon)
					codon_l[index] = "C"
					new_codon = ''.join(codon_l)
					if not new_codon in stop_codons:
						if DNA_CODON_TABLE[new_codon] == aa:
							syn_nonsyn.append("s")
						else:
							syn_nonsyn.append("ns")
			syn_sites = 3*(syn_nonsyn.count("s")/len(syn_nonsyn))
			nonsyn_sites = 3*(syn_nonsyn.count("ns")/len(syn_nonsyn))
			syn_sitesl.append(syn_sites)
			nonsyn_sitesl.append(nonsyn_sites)
	return sum(syn_sitesl), sum(nonsyn_sitesl)	
#*******************************************************************************	 
# Return count of fourfold SNPs	 
#*******************************************************************************
def fourfold_SNPs(codon):		
	fourfold = 0
	amino_acids = []
	if not "N" in codon:
		for base in codon:
			if base in IUPAC_code.keys():
				SNP_index = codon.index(base)
				for nuc in nucs:
					new_codon = codon.replace(codon[SNP_index], nuc)
					amino_acids.append(DNA_CODON_TABLE[new_codon])
		if len(list(set(amino_acids))) == 1:
			if amino_acids[0] in four_fold_codons.keys():
				fourfold += 1
	return fourfold						
#*******************************************************************************
# Return count of zerofold SNPs
#*******************************************************************************
def zerofold_SNPs(codon):	
	zerofold = 0
	amino_acids = []
	if not "N" in codon:
		for base in codon:
			if base in IUPAC_code.keys():
				SNP_index = codon.index(base)
				for nuc in nucs:
					new_codon = codon.replace(codon[SNP_index], nuc)
					amino_acids.append(DNA_CODON_TABLE[new_codon])
		if len(list(set(amino_acids))) == 4:
			zerofold += 1
	return zerofold
#*******************************************************************************
# Return count of SW fourfold SNPs							
#*******************************************************************************
def fourfold_SW(codon):		
	fourfold = 0
	amino_acids = []
	if not "N" in codon:
		for base in codon:
			if base == "S" or base == "W":
				SNP_index = codon.index(base)
				for nuc in nucs:
					new_codon = codon.replace(codon[SNP_index], nuc)
					amino_acids.append(DNA_CODON_TABLE[new_codon])
		if len(list(set(amino_acids))) == 1:
			if amino_acids[0] in four_fold_codons.keys():
				fourfold += 1
	return fourfold
#*******************************************************************************	
# Return count of SW zerofold SNPs	
#*******************************************************************************
def zerofold_SW(codon):	
	zerofold = 0
	amino_acids = []
	if not "N" in codon:
		for base  in codon:
			if base == "S" or base == "W":
				SNP_index = codon.index(base)
				for nuc in nucs:
					new_codon = codon.replace(codon[SNP_index], nuc)
					amino_acids.append(DNA_CODON_TABLE[new_codon])
		if len(list(set(amino_acids))) == 4:
			zerofold += 1
	return zerofold		        
#*******************************************************************************
# Test functions
#*******************************************************************************
def test_count_sites():
	if not count_sites('ATG')[1] == 3:
		raise Exception
	# Add other arguments here for test
test_count_sites()

def test_fourfold():
	if not count_fourfold_sites('TYT') == 0.5:
		raise Exception
	elif not count_fourfold_sites('TTY') == 0:
		raise Exception
	elif not count_fourfold_sites('SGG') == 1:
		raise Exception
test_fourfold()

def test_zerofold():
	if not count_zerofold_sites('ATG') == 3:
		raise Exception
	elif not count_zerofold_sites('TTY') == 2:
		raise Exception
	elif not count_zerofold_sites('TGA') == 0:
		raise Exception
test_zerofold()

