#!/usr/bin/python
import sys

for pop in "CZC", "CZP", "E", "H", "I", "OC", "OP", "SP":
	for i in range(1, 28):
		print "#!/bin/bash -l","\n",\
		"#SBATCH --mail-user homa.papoli@ebc.uu.se","\n",\
		"#SBATCH -A b2010010","\n",\
		"#SBATCH -p core","\n",\
		"#SBATCH -n 1","\n",\
		"#SBATCH -t 24:00:00","\n",\
		"#SBATCH -J Rho","\n",\

		print "./SNP_Rho_BalSel_V3.py pop_snpPair/${chr}_col.${pop}.snpPair_runALL.txt $snpfile > ${pop}/balletFreqInput_${pop}_15Oct2015_${chr}.rho"