#!/usr/bin/python3 

import numpy as np
import sys

n_pop = 190

file_in = sys.argv[1] # *.jm
gt_array = np.zeros(n_pop)
minN = n_pop
repGT = np.zeros(n_pop)
strPreChr = ''
Outfile = open(file_in+'repgt.txt','w')
for line in open(file_in):

	cell = line.strip().split('\t')
	strPos = cell[0]
	strChr = '_'.join(strPos.split('_')[0:2])

	gt_array = np.asarray(cell[1:])
	print(strPos,repGT,len(repGT))
	print(strPos,gt_array,len(gt_array))
	deter_array = np.array((repGT == gt_array),dtype=int)
	if strPreChr == '':
		repGT = gt_array
		repPos = strPos
		minN = n_pop
		pass
	elif strChr == strPreChr:
		if sum(deter_array) / n_pop > 0.95:
			if (gt_array == 'n').sum() < minN:
				minN = (gt_array == 'n').sum()
				repPos = strPos 
				repGT = gt_array
			else:
				pass
		else:
			print(repPos, '\t'.join(repGT),sep='\t',file=Outfile)
			repPos = strPos
			repGT = gt_array 
			minN = n_pop
	else:
		print(repPos, '\t'.join(repGT),sep='\t',file=Outfile)
		repPos = strPos
		repGT = gt_array
		minN = n_pop
	strPrePos = cell[0]
	strPreChr = '_'.join(strPos.split('_')[0:2])	
	
