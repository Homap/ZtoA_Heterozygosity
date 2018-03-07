#!/usr/bin/python
import sys
from fasta import readfasta

# To get the count of transitions and transversions as a sanity check of SNP data
# Run: ./count_trans_trans.py fasta.fa > trans_transv.txt

def trans_tranv_count(sequence):
	"""Return heterozygosity from
		any sequence as a string """
	# The first thing is to test the input
	if not isinstance(sequence, str):
		raise Exception("Sequence is not a string")
	R = len([base.upper() for base in sequence if base.upper()=="R"])
	Y = len([base.upper() for base in sequence if base.upper()=="Y"])
	S = len([base.upper() for base in sequence if base.upper()=="S"])
	W = len([base.upper() for base in sequence if base.upper()=="W"])
	K = len([base.upper() for base in sequence if base.upper()=="K"])
	M = len([base.upper() for base in sequence if base.upper()=="M"])
	
	return (R, Y, S, W, K, M) 
	
################################################################################
fastafile = open(sys.argv[1], "r")

fastadict = readfasta(fastafile)

R_l = [] 
Y_l = [] 
S_l = [] 
W_l = [] 
K_l = [] 
M_l = []
for key in fastadict.keys():
	R = trans_tranv_count(fastadict[key])[0]
	Y = trans_tranv_count(fastadict[key])[1]
	S = trans_tranv_count(fastadict[key])[2]
	W = trans_tranv_count(fastadict[key])[3]
	K = trans_tranv_count(fastadict[key])[4]
	M = trans_tranv_count(fastadict[key])[5]
	R_l.append(R)
	Y_l.append(Y)
	S_l.append(S)
	W_l.append(W)
	K_l.append(K)
	M_l.append(M)
	
#print ("R"+"\t"+"Y"+"\t"+"S"+"\t"+"W"+"\t"+"K"+"\t"+"M"+"\t"+"Transition"+"\t"+"Transversion") 
print (str(sum(R_l))+"\t"+str(sum(Y_l))+"\t"+str(sum(S_l))+"\t"+str(sum(W_l))+"\t"+str(sum(K_l))+"\t"+str(sum(M_l))+"\t"+str(sum(R_l)+sum(Y_l))+"\t"+str(sum(S_l)+sum(W_l)+sum(K_l)+sum(M_l)))

	
