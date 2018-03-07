#!/usr/bin/python
import sys

"""
Written by Homa Papoli, Nov 2017
scaf_chr_het.py takes a chromosome, scaffold file and a heterozygosity file. 
If heterozygosity is from CDS, keyword CDS must be provided as the third 
argument, otherwise, intron or intergenic.
File structure:
Chromosome file: chromosome"\t"scaffold
Heterozygosity file for intron and intergenic:
scaffold	start	end	SNPs	bases	heterozygosity
Heterozygosity file for CDS:
scaffold_geneID	Pi_syn	syn_sites	Pi_syn/syn_sites	Pi_nonsyn	nonsyn_sites	Pi_nonsyn/nonsyn_sites	length_orf
USAGE:
./scaf_chr_het.py all_chr Apaloderma_vittatum.all.intron.het intron
"""


chrfile = open(sys.argv[1], "r")
hetfile = open(sys.argv[2], "r")
functional = sys.argv[3]

if functional == "CDS":
	# Read chrfile into a dictionary
	chrdict = {}
	for line in chrfile:
		line = line.strip("\n").split("\t")
		key, value = line[0], line[1]
		if key in chrdict.keys():
			chrdict[key].append(value)
		else:
			chrdict[key] = [value]
	# Read het file into a dictionary
	hetdict = {}
	for line in hetfile:
		line = line.strip("\n").split("\t")
		if not line[0].startswith("geneID"):
			key, value = line[0].split("_")[0], [line[1], line[2], line[3], line[4], line[5], line[6], line[7], line[8]]
			if key in hetdict.keys():
				hetdict[key].append(value)
			else:
				hetdict[key] = [value]
#	print hetdict
#	print "\t".join(["CHROM", "Pi_syn", "syn_sites", "Pi_nonsyn", "nonsyn_sites", "fourfold", "four_fold_sites", "zerofold", "zero_fold_sites"])
	for chrom in chrdict.keys():
		Pi_syn = 0
		syn_sites = 0
		Pi_nonsyn = 0
		nonsyn_sites = 0
		fourfold = 0
		four_fold_sites = 0
		zerofold = 0
		zero_fold_sites = 0
		for scaffold in chrdict[chrom]:
			if scaffold in hetdict.keys():
				for coord in hetdict[scaffold]:
					Pi_syn += float(coord[0])
					syn_sites += float(coord[1])
					Pi_nonsyn += float(coord[2])
					nonsyn_sites += float(coord[3])
					fourfold += float(coord[4])
					four_fold_sites += float(coord[5])
					zerofold += float(coord[6])
					zero_fold_sites += float(coord[7])
		print "\t".join([chrom, str(Pi_syn), str(syn_sites), str(Pi_nonsyn), str(nonsyn_sites), str(fourfold), str(four_fold_sites), str(zerofold), str(zero_fold_sites)])
else:
	# Read chrfile into a dictionary
	chrdict = {}
	for line in chrfile:
		line = line.strip("\n").split("\t")
		key, value = line[0], line[1]
		if key in chrdict.keys():
			chrdict[key].append(value)
		else:
			chrdict[key] = [value]
	#print chrdict
	# Read het file into a dictionary
	hetdict = {}
	for line in hetfile:
		line = line.strip("\n").split("\t")
		key, value = line[0].split("_")[0], [line[3], line[4]]
		if key in hetdict.keys():
			hetdict[key].append(value)
		else:
			hetdict[key] = [value]
	#print hetdict
#	print "\t".join(["CHROM", "SNP", "TOTAL"])
	for chrom in chrdict.keys():
		SNP = 0
		total = 0
		for scaffold in chrdict[chrom]:
			if scaffold in hetdict.keys():
				for coord in hetdict[scaffold]:
					if not coord[0] == "NA":
						SNP = SNP + int(coord[0])
					if not coord[1] == "NA":
						total = total + int(coord[1])
		print "\t".join([chrom, str(SNP), str(total)])
	

	
	
