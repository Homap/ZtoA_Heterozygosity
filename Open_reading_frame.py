#!/usr/bin/python
from __future__ import division
import sys
import collections
from fasta import readfasta
from SynNonSynCounts import count_sites
from SynNonSynCounts import count_fourfold_sites
from SynNonSynCounts import count_zerofold_sites

##########################################################################################
# Written by Homa Papoli - 16 May 2016 
# Updated Sept 2017                                              
# This script contains the following functions:
# revcomp(dna, reverse=True, complement=True) : Produce a reverse complement strand of DNA
# dna_translate(cdna, code = DNA_CODON_TABLE) : Produce a protein sequence from DNA
# findlongest(sequencelist, fragonly=True) : Find the longest protein sequence
# Usage: 
# If run: ./Open_reading_frame.py CDS.fa all > output
# Script counts all synonymous and nonsynonymous polymorphisms
# If run: ./Open_reading_frame.py CDS.fa SW > output
# The script counts only S and W sites. 
# CDS.fa is produced by extract_CDS_from_fasta.py and contains a fasta sequence with IUPAC
# coded SNPs.
##########################################################################################

#*******************
# Specify the inputs
#*******************
inFile = open(sys.argv[1], 'r')
argument = sys.argv[2]

#************************************************************
# Read the fasta file into a dictionary
#************************************************************
fastaseq = readfasta(inFile)
# Specify four nucleotides
#************************************************************
nucs = ["A", "T", "C", "G"]
#************************************************************
# Specify IUPAC codes as a dictionary
#************************************************************
IUPAC_code = {'R':['A', 'G'], 'Y':['C', 'T'], 'S':['G', 'C'], 'W':['A', 'T'], 'K':['G', 'T'], 'M':['A', 'C']}#, 'B':['C', 'G', 'T'], 'D':['A', 'G', 'T'], 'H':['A', 'C', 'T'], 'V':['A', 'C', 'G']}
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
    
#*************************************************************
# Define a function to reverse compliment a strand
#*************************************************************
def revcomp(dna, reverse=True, complement=True):
	""" Takes a sequence of DNA and converts it to its
	reverse compliment, if only compliment is intended
	set reverse to False """
	bases = 'ATGCatgcWSRYMKwsrymkHBVDhbvdNnTACGTACGWSYRKMWSYRKMDVBHDVBHNN'
	complement_dict = {} # build a dictionary that contains each base with its complement.
	for i in range(30):
		complement_dict[bases[i]] = bases[i+30]
	if reverse: # if reverse is True, default is true
		dna = reversed(dna)
		result_as_list = None # define an empty list
	if complement: # if complement is true, default is true
		result_as_list = [complement_dict[base] for base in dna]
	else:
		result_as_list = [base for base in dna]
	return ''.join(result_as_list)
	
# Test revcomp() function for two cases:	
def test_revcomp():
	if not revcomp("ATC") == "GAT":
		raise Exception
	elif not revcomp("ARTTT") == "AAAYT":
		raise Exception

test_revcomp()
	
#************************************************************** 
# Define a function to translate DNA into 6 reading frames
#**************************************************************	
def dna_translate(cdna, code = DNA_CODON_TABLE):
	""" Takes a CDS and translate it in the 6 possible
	open reading frames. 3 frames 5'->3' and 3 frames
	3'->5'. It returns a tuple, containing two lists. 
	The first is a list of all possible proteins and the
	second is a list of all possible open reading frames."""
	allprot = []
	dna_orf = []
	# X is returned when an IUPAC code corresponds to a nonsynonymous change.
	allprot.append("".join([ code.get(cdna[i:i+3], "X") for i in xrange(0, len(cdna), 3) ]) )
	dna_orf.append("".join([ cdna[i:i+3] for i in xrange(0, len(cdna), 3) ]) )
	allprot.append("".join([ code.get(cdna[i:i+3], "X") for i in xrange(1, len(cdna), 3) ]) )
	dna_orf.append("".join([ cdna[i:i+3] for i in xrange(1, len(cdna), 3) ]) )
	allprot.append("".join([ code.get(cdna[i:i+3], "X") for i in xrange(2, len(cdna),3) ]) )
	dna_orf.append("".join([ cdna[i:i+3] for i in xrange(2, len(cdna), 3) ]) )
	revcdna = revcomp(cdna)
	allprot.append("".join([ code.get(revcdna[i:i+3], "X") for i in xrange(0,len(revcdna),3) ]) )
	dna_orf.append("".join([ revcdna[i:i+3] for i in xrange(0, len(revcdna), 3) ]) )
	allprot.append("".join([ code.get(revcdna[i:i+3], "X") for i in xrange(1,len(revcdna),3) ]) )
	dna_orf.append("".join([ revcdna[i:i+3] for i in xrange(1, len(revcdna), 3) ]) )
	allprot.append("".join([ code.get(revcdna[i:i+3], "X") for i in xrange(2,len(revcdna),3) ]) )
	dna_orf.append("".join([ revcdna[i:i+3] for i in xrange(2, len(revcdna), 3) ]) )
	return allprot, dna_orf

#*************************************************************************
# Define a function to find the longest continuous ORF 
# I want the script to give me the DNA sequence of the longest ORF as well
#*************************************************************************
def findlongest(sequencelist, fragonly=True):
	""" This function finds the longest open reading frame
	and return that."""
	# finds the longest continuous ORF
	lenscores={}
	fragscores=[]
	bestorflist=[]
	bestorfdna=[]
	for index, seq in enumerate(sequencelist[0]):
		orfs=seq.rstrip().split('_')
		print orfs
		lenlist=[len(frag) for frag in orfs]
		print lenlist
		lenscores[index] = (max(lenlist))
		print lenscores	
		if fragonly:
			bestorflist.append( orfs[lenlist.index(max(lenlist))] )
			print bestorflist
		fragscores.append(len(lenlist)-lenlist.count(0)) # lenlist.count(0) because if there are two adjacent stop codons, it becomes '' after split("_"). e.g: QAVVFLFLLPPRNL__IHRKHQLSQGEICISPKHLHRSYLD
		print len(lenlist), lenlist.count(0), fragscores
	maxlenlist = [max(lenscores, key=lambda i:lenscores[i])] # get the sequence with the longest ORF
	masterlist = list(set(maxlenlist))	
	
	if (len(masterlist)==1):
		if fragonly: # fragonly is by defaul True. When True, findlongest() returns both protein and DNA sequence of the reading frame. If set to False, it only returns the protein sequence.
			return [ [bestorflist[ masterlist[0] ]], sequencelist[1][masterlist[0]], masterlist[0] ]
		else:
			return [ [sequencelist[0][ masterlist[0] ]], masterlist[0] ]
	
	else:
		if (max(fragscores) > 1): # this means there are more than one segment when cDNA is translated, stopped by multiple stop codons. 
			if fragonly:
				return [[bestorflist[x] for x in masterlist],  masterlist]
			else:
				return [[sequencelist[0][x] for x in masterlist], masterlist]
		else:
			return [["---"], [0]]		
			
#***************************************************************************
# Run the functions for the CDS fasta sequence
#***************************************************************************
# Define dna_orf dictionary
dna_orf = {}
for key, value in fastaseq.iteritems():
#	print key, value
	protein = dna_translate(value)
	print protein
	orfdna = findlongest(protein)
	print orfdna
	dna_orf[key] = orfdna[1]
print " Open reading frame", dna_orf

## Go through each sequence in the dna_orf dictionary 
## and add the codons with IUPAC code in the list for each gene.
#IUPAC_aa = {} # IUPAC_aa is a dictionary of codons containing an IUPAC coded base
#SNs_counts = {}
#for gene, sequence in dna_orf.iteritems():
##	print gene, sequence
#	SNs_counts[gene] = list(count_sites(sequence))
##	print SNs_counts
#	sequence = sequence.upper()
#	last_codon_start = len(sequence) - 2
#	for base in range(0, last_codon_start, 3):
#		codon = sequence[base:base+3]
#		IUPAC_set = IUPAC_code.keys()
#		matching = [i for i in codon if i in IUPAC_set] # check if codon contains a SNP
#		if len(matching) < 2:
#			if gene in IUPAC_aa.keys():
#				IUPAC_aa[gene].append(codon)
#			else:
#				IUPAC_aa[gene] = [codon]
#						
##print "IUPAC_aa", IUPAC_aa
## Here I count the number of synonymous and nonsynonymous polymorphisms
## If the polymorphism changes the amino acid, it is synonymous and if not
## it is nonsynonymous.
#if argument == "all":
#	header = ["geneID", "Pi_syn", "syn_sites", "Pi_syn/syn_sites", "Pi_nonsyn", "nonsyn_sites", "Pi_nonsyn/nonsyn_sites", "fourfold", "four_fold_sites", "fourfold/four_fold_sites", "zerofold", "zero_fold_sites", "zerofold/zero_fold_sites", "length_orf"]
#	print '\t'.join(header)	
#	for key, value in IUPAC_aa.iteritems():
#		syn_sites = SNs_counts[key][0] # synonymous sites
#		Nonsyn_sites = SNs_counts[key][1] # non synonymous sites
#		four_fold_sites = SNs_counts[key][2] # four fold degenerate sites
#		zero_fold_sites = SNs_counts[key][3] # zero fold degenerate sites
##		print key, value, syn_sites, four_fold_sites
#		synonymous = 0 # synonymous SNPs
#		nonsynonymous = 0 # nonsynonymous SNPs
#		fourfold = 0 # fourfold SNPs
#		zerofold = 0 # zerofold SNPs
#		for codon in value:
#			for base in codon:
#				if base in IUPAC_code.keys():
##					print base
#					# Get index of the variable site
#					SNP_index = codon.index(base)
#					# Define a list to store all 4 amino acids
#					# Define a loop to replace SNP with all four nucleotides
##					zero_fold_evaluation = [DNA_CODON_TABLE[codon.replace(codon[SNP_index], i)] for i in nucs]
##					if len(list(set(zero_fold_evaluation))) == 4:
##						zerofold += 1
##					print zero_fold_evaluation
##					print key, zerofold
#					codon1 = codon.replace(base, IUPAC_code[base][0])
#					codon2 = codon.replace(base, IUPAC_code[base][1])
#					aa1 = DNA_CODON_TABLE.get(codon1,"X") 
#					aa2 = DNA_CODON_TABLE.get(codon2, "X")
##					print codon1, codon2, aa1, aa2
#					if aa1 == aa2:
#						if aa1 != "X":
#							if aa1 in four_fold_codons.keys():
#								synonymous += 1
#								#***************************************************
#								fourfold += 1 
##								print codon, codon1, codon2, aa1, aa2, 'synonymous', synonymous, fourfold
#							elif aa1 not in four_fold_codons.keys():
#								synonymous += 1
##								print codon, codon1, codon2, aa1, aa2, 'synonymous', synonymous, fourfold
#					else:
#						nonsynonymous += 1
##						print codon, codon1, codon2, aa1, aa2, 'nonsynonymous', nonsynonymous, fourfold
##		print syn_sites				
#		if syn_sites == 0 and Nonsyn_sites != 0:
#			out = [key,synonymous,syn_sites,"NA",nonsynonymous,Nonsyn_sites,nonsynonymous/Nonsyn_sites, "NA", "NA", "NA", len(dna_orf[key])]
#			print '\t'.join(map(str,out))
#		elif syn_sites != 0 and Nonsyn_sites == 0:
#			if four_fold_sites == 0:
#				out = [key,synonymous,syn_sites,synonymous/syn_sites,nonsynonymous,Nonsyn_sites,"NA", "NA", "NA", "NA", len(dna_orf[key])]
#				print '\t'.join(map(str,out))
#			else:
#				out = [key,synonymous,syn_sites,synonymous/syn_sites,nonsynonymous,Nonsyn_sites,"NA", fourfold, four_fold_sites, fourfold/four_fold_sites, len(dna_orf[key])]
#				print '\t'.join(map(str,out))
#		elif syn_sites == 0 and Nonsyn_sites == 0:
#			out = [key,synonymous,syn_sites,"NA",nonsynonymous,Nonsyn_sites,"NA", "NA", "NA", "NA", len(dna_orf[key])]
#			print '\t'.join(map(str,out))
#		else:
#			if four_fold_sites == 0:
#				out = [key,synonymous,syn_sites,synonymous/syn_sites,nonsynonymous,Nonsyn_sites,"NA", "NA", "NA", "NA", len(dna_orf[key])]
#				print '\t'.join(map(str,out))
#			else:
#				out = [key,synonymous,syn_sites,synonymous/syn_sites,nonsynonymous,Nonsyn_sites,"NA", fourfold, four_fold_sites, fourfold/four_fold_sites, len(dna_orf[key])]
#				print '\t'.join(map(str,out))
#else:
#	header = ["geneID", "Pi_syn", "syn_sites", "Pi_syn/syn_sites", "Pi_nonsyn", "nonsyn_sites", "Pi_nonsyn/nonsyn_sites", "fourfold", "four_fold_sites", "fourfold/four_fold_sites", "zerofold", "zero_fold_sites", "zerofold/zero_fold_sites", "length_orf"]
#	print '\t'.join(header)		
#	for key, value in IUPAC_aa.iteritems():
#		syn_sites = SNs_counts[key][0]
#		Nonsyn_sites = SNs_counts[key][1]
#		four_fold_sites = SNs_counts[key][2] # four fold degenerate sites
#		zero_fold_sites = SNs_counts[key][3] # zero fold degenerate sites
##		print key, value, syn_sites
#		synonymous = 0
#		nonsynonymous = 0
#		fourfold = 0 # fourfold SNPs
#		zerofold = 0 # zerofold SNPs
#		for codon in value:
#			for base in codon:
#				if base in IUPAC_code.keys():
#					if base == "W" or base == "S":
#					# Get index of the variable site
#						SNP_index = codon.index(base)
#						# Define a list to store all 4 amino acids
#						# Define a loop to replace SNP with all four nucleotides
#						zero_fold_evaluation = [DNA_CODON_TABLE[codon.replace(codon[SNP_index], i)] for i in nucs]
#						if len(list(set(zero_fold_evaluation))) == 4:
#							zerofold += 1
#						print zero_fold_evaluation
#						print key, zerofold
##						print codon, base
#						codon1= codon.replace(base, IUPAC_code[base][0])
#						codon2= codon.replace(base, IUPAC_code[base][1])
#						aa1 = DNA_CODON_TABLE.get(codon1,"X") 
#						aa2 = DNA_CODON_TABLE.get(codon2, "X")
#						if aa1 == aa2:
#							if aa1 in four_fold_codons.keys():
#								synonymous += 1
#								#***************************************************
#								fourfold += 1 
##								print codon, codon1, codon2, aa1, aa2, 'synonymous', synonymous, fourfold
#							elif aa1 not in four_fold_codons.keys():
#								synonymous += 1
##								print codon, codon1, codon2, aa1, aa2, 'synonymous', synonymous, fourfold
#						else:
#							nonsynonymous += 1
##							print codon, codon1, codon2, aa1, aa2, 'nonsynonymous', nonsynonymous, fourfold
##			print syn_sites				
#			if syn_sites == 0 and Nonsyn_sites != 0:
#				out = [key,synonymous,syn_sites,"NA",nonsynonymous,Nonsyn_sites,nonsynonymous/Nonsyn_sites, "NA", "NA", "NA", len(dna_orf[key])]
#				print '\t'.join(map(str,out))
#			elif syn_sites != 0 and Nonsyn_sites == 0:
#				if four_fold_sites == 0:
#					out = [key,synonymous,syn_sites,synonymous/syn_sites,nonsynonymous,Nonsyn_sites,"NA", "NA", "NA", "NA", len(dna_orf[key])]
#					print '\t'.join(map(str,out))
#				else:
#					out = [key,synonymous,syn_sites,synonymous/syn_sites,nonsynonymous,Nonsyn_sites,"NA", fourfold, four_fold_sites, fourfold/four_fold_sites, len(dna_orf[key])]
#					print '\t'.join(map(str,out))
#			elif syn_sites == 0 and Nonsyn_sites == 0:
#				out = [key,synonymous,syn_sites,"NA",nonsynonymous,Nonsyn_sites,"NA", "NA", "NA", "NA", len(dna_orf[key])]
#				print '\t'.join(map(str,out))
#			else:
#				if four_fold_sites == 0:
#					out = [key,synonymous,syn_sites,synonymous/syn_sites,nonsynonymous,Nonsyn_sites,"NA", "NA", "NA", "NA", len(dna_orf[key])]
#					print '\t'.join(map(str,out))
#				else:
#					out = [key,synonymous,syn_sites,synonymous/syn_sites,nonsynonymous,Nonsyn_sites,"NA", fourfold, four_fold_sites, fourfold/four_fold_sites, len(dna_orf[key])]
#					print '\t'.join(map(str,out))
	 
