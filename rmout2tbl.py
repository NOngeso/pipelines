#!/usr/bin/python3

import sys


dicRN2Annot = {}
file_soyte = 'redbean.repeatlib.fa.bn.1e10.SoyBase_TE_Fasta.txt.parsed.fa.out'
for line in open(file_soyte):
	if line[0] == '#':
		continue
	cell  = line.strip().split('\t')
	if len(cell[0].split('#')) < 2 or len(cell[1].split('#')) < 2:
		continue
	strRN = cell[0].split('#')[0]
	strTE = cell[1].split('#')[1]
	dicRN2Annot[strRN] = strTE.strip()


file_rmout = 'adzuki.ref.fa.out'
Outfile    = open(file_rmout+'.table_ver2','w')
Outfile_Rin = open(file_rmout+'.Rin','w')
print('Chromosome','Position','TE',sep='\t',file=Outfile_Rin)
dicRC2lengSum = {}
i = 0
for line in open(file_rmout):
	if i < 3 :
		i += 1
		continue
	cell = line.strip().split()
	strC = cell[4]
	strP = cell[5]
	strRN  = cell[9]
	strRC  = cell[10] # Repeat class
	try:
		if 'caf' not in strC:
			print(strC, strP, dicRN2Annot[strRN],sep='\t',file=Outfile_Rin)
	except KeyError:
		pass
	try:
		strRC	= dicRN2Annot[strRN]
	except KeyError:
		strRC	= strRC
	#if 'LTR/Others' in strRC:
	#	try:
	#		strRC_alt = dicRN2Annot[strRN]
	#	except KeyError:
	#		strRC_alt = ''
	#	if 'LTR' in strRC_alt:
	#		strRC = strRC_alt
	strTP  = cell[8]
	if strTP == '+':
		intLen = abs(int(cell[12].strip()) - int(cell[11].strip())) + 1 
	elif strTP == 'C':
		intLen = abs(int(cell[13].strip()) - int(cell[12].strip())) + 1
	try:
		dicRC2lengSum[strRC][0] += intLen
		dicRC2lengSum[strRC][1] += 1
	except KeyError:
		dicRC2lengSum[strRC] = [intLen,1]
for strRC in dicRC2lengSum:
	print(strRC,dicRC2lengSum[strRC][0],dicRC2lengSum[strRC][1],sep='\t',file=Outfile)
Outfile.close()
