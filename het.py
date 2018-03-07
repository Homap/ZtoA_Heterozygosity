#!/usr/bin/python
from __future__ import division
import sys

def heterozygosity(sequence):
	"""Return heterozygosity from
		any sequence as a string """
	# The first thing is to test the input
	if not isinstance(sequence, str):
		raise Exception("Sequence is not a string")
	homozygote='ATCG'
	SNPs='RYSWKMBDHV'	
	homozygote_count = len([base.upper() for base in sequence if base.upper() in homozygote])		
	SNPs_count = len([base.upper() for base in sequence if base.upper() in SNPs])		
	if not homozygote_count+SNPs_count == 0:
		return SNPs_count, homozygote_count+SNPs_count, (SNPs_count / (homozygote_count+SNPs_count))
	else:
		return 'NA', 'NA', 'NA'


def test():
	#check heterzygosity
	if not heterozygosity('ATWK')[2] == 0.5:
		raise Exception
test()
	
	
