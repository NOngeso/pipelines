#!/usr/bin/python3
import sys
for line in sys.stdin:
	cell = line.strip().split('\t')
	gt_list = cell[3].split(',')
	if gt_list.count('N') <= 20:
		pass
	else : continue
	dicGT2num = {}
	for gt in gt_list:
		if gt == 'N':
			continue
		try:
			dicGT2num[gt] += 1
		except KeyError:
			dicGT2num[gt] = 1
	total = sum(list(dicGT2num.values()))
	dicGT2num_list = list(dicGT2num.keys())
	dicGT2num_list.sort(key=lambda x : dicGT2num[x])	
	if dicGT2num[dicGT2num_list[0]]/total < 0.1:# minor allele freq
		continue
	if 'N' in gt_list:
		if len(set(gt_list)) >  2:
			print(line.strip())
	else: 
		if len(set(gt_list)) > 1:
			print(line.strip())
