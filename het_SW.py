#!/usr/bin/python
from __future__ import division
import sys

def heterozygositySW(sequence):
	"""Return heterozygosity from
		any sequence as a string """
	# The first thing is to test the input
	if not isinstance(sequence, str):
		raise Exception("Sequence is not a string")
	homozygote='ATCG'
	SNPs='RYSWKMBDHV'
	SWs='SW'	
	homozygote_count = len([base.upper() for base in sequence if base.upper() in homozygote])		
	SNPs_count = len([base.upper() for base in sequence if base.upper() in SNPs])
	SWs_count = len([base.upper() for base in sequence if base.upper() in SWs])		
	if not homozygote_count+SNPs_count == 0:
		return SWs_count, homozygote_count+SNPs_count, (SWs_count / (homozygote_count+SNPs_count))
	else:
		return 'NA', 'NA', 'NA'


def test():
	#check heterzygosity
	if not heterozygositySW('ATGSWS')[2] == 0.5:
		raise Exception
test()
