#!/usr/bin/python
from __future__ import division
from datetime import datetime
import sys

fasta = open(sys.argv[1], "r")
#heterozygosity = open(sys.argv[2], "w")
def heterozygosity(fasta):
	total_A = 0
	total_T = 0
	total_C = 0
	total_G = 0
	total_R = 0
	total_Y = 0
	total_S = 0
	total_W = 0
	total_K = 0
	total_M = 0

	for line in fasta:
		if line.startswith(">"):
			print line.rstrip('\n')[1:]
		else:
			line = line.rstrip('\n').upper()
			A = line.count('A')
			T = line.count('T')
			C = line.count('C')
			G = line.count('G')
			R = line.count('R')
			Y = line.count('Y')
			S = line.count('S')
			W = line.count('W')
			K = line.count('K')
			M = line.count('M')
			total_A = total_A + A
			total_T = total_T + T
			total_C = total_C + C
			total_G = total_G + G
			total_R = total_R + R
			total_Y = total_Y + Y
			total_S = total_S + S
			total_W = total_W + W
			total_K = total_K + K
			total_M = total_M + M
	
		heterozygosity.All_bases = total_A + total_T + total_C + total_G + total_R + total_Y + total_S + total_W + total_K + total_M 
		heterozygosity.Heterozygote_sites = total_R + total_Y + total_S + total_W + total_K + total_M
		if heterozygosity.All_bases!= 0:
			return heterozygosity.Heterozygote_sites / heterozygosity.All_bases
	
#print "Heterozygosity: " + str(heterozygosity(fasta))
#print "Total number of bases: " + str(heterozygosity.All_bases/1000000)
#print "Number of heterozygote sites: " + str(Heterozygote_sites)
print str(round(heterozygosity(fasta), 4)) + "\t" + str(heterozygosity.All_bases/1000000) + "\t" + str(heterozygosity.Heterozygote_sites)

fasta.close()
