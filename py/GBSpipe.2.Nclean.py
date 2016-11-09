#!/usr/bin/python3
import sys
for line in sys.stdin:
	cell = line.strip().split('\t')
	gt_list = cell[3].split(',')
	if gt_list.count('N') <= 20:
		pass
	else : continue
	
	if 'N' in gt_list:
		if len(set(gt_list)) > 2:
			print(line.strip())
	else: 
		if len(set(gt_list)) > 1:
			print(line.strip())
