#!/usr/bin/python
from __future__ import division
import sys
import collections
from fasta import readfasta
from SynNonSynCounts import count_sites

#*******************
# Specify the inputs
#*******************
inFile = open(sys.argv[1], 'r')

#************************************************************
# Read the fasta file into a dictionary
#************************************************************
fastaseq = readfasta(inFile)
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

# Do not allow for mutations into stop codons (Yang)
stop_codons = ['TGA', 'TAA', 'TAG', 'TRA', 'TAR']
def count_sites(sequence, nucs=nucs, DNA_CODON_TABLE=DNA_CODON_TABLE):
	syn_sitesl = [] # list of synonymous sites
	nonsyn_sitesl = [] # list of non-synonymous sites
	four_fold = [] # list of four-fold degenerate sites
	for i in range(0, len(sequence)-2, 3): # a codon every three bases
		syn_nonsyn = []
		codon = sequence[i:i+3]
		if not 'N' in codon:
			codonl = list(codon) # convert codon string into list
			IUPAC_set = IUPAC_code.keys()
			matching = [i for i in codonl if i in IUPAC_set] # check if codon contains a SNP
			if matching: # if there is a SNP in codon
				if len(matching) == 1: # if there is only one SNP in a codon
					syn_nonsyn1 = []
					syn_nonsyn2 = []
					# *************************
					SNP = list(matching)[0]
					codon1 = codon.replace(SNP, IUPAC_code[SNP][0])
					aa1 = DNA_CODON_TABLE[codon1]
					codon2 = codon.replace(SNP, IUPAC_code[SNP][1])
					aa2 = DNA_CODON_TABLE[codon2]
					# Synonymous-Nonsynonymous sites
					for index in range(0,3,1):
						for nuc in nucs:
							if not codon1[index] == nuc:
								codon1_l = list(codon1)
								codon1_l[index] = nuc
								new_codon1 = ''.join(codon1_l)
								if not new_codon1 in stop_codons:
									if DNA_CODON_TABLE[new_codon1] == aa1:
										syn_nonsyn1.append("s")
									else:
										syn_nonsyn1.append("ns")
							if not codon2[index] == nuc:
								codon2_l = list(codon2)
								codon2_l[index] = nuc
								new_codon2 = ''.join(codon2_l)
								if not new_codon2 in stop_codons:
									if DNA_CODON_TABLE[new_codon2] == aa2:
										syn_nonsyn2.append("s")
									else:
										syn_nonsyn2.append("ns")
					syn_sites1 = 3*(syn_nonsyn1.count("s")/len(syn_nonsyn1))
					nonsyn_sites1 = 3*(syn_nonsyn1.count("ns")/len(syn_nonsyn1))
					syn_sites2 = 3*(syn_nonsyn2.count("s")/len(syn_nonsyn2))
					nonsyn_sites2 = 3*(syn_nonsyn2.count("ns")/len(syn_nonsyn2))
					SNP_syn_sites = (syn_sites1+syn_sites2)/2
					SNP_nonsyn_sites = (nonsyn_sites1+nonsyn_sites2)/2
					syn_sitesl.append(SNP_syn_sites)
					nonsyn_sitesl.append(SNP_nonsyn_sites)
					# 4-fold degenerate sites with SNP *************************
					if aa1 == aa2:
						if aa1 in four_fold_codons.keys():
							four_fold.append(1)
					else:
						if aa1 in four_fold_codons.keys() or aa2 in four_fold_codons.keys():
							four_fold.append(0.5)
					# **********************************************************
			else:
				aa = DNA_CODON_TABLE[codon]
				for index in range(0,3,1):					
					for nuc in nucs:
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
				aa = DNA_CODON_TABLE[codon]
				# 4-fold degenerate sites no SNP ****************************
				if DNA_CODON_TABLE[codon] in four_fold_codons.keys():
					four_fold.append(1)
						
	return sum(syn_sitesl), sum(nonsyn_sitesl), sum(four_fold)
	         
	         
seq = "ATGCCCCCYCYC"
print count_sites(seq)
