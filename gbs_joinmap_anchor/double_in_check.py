#!/usr/bin/python3

import sys
done = []
file_jo = sys.argv[1]#'joinmap_ordered.txt.sorted.curated'
for line in open(file_jo):
	if line[0] == '#':
		continue
	cell = line.strip().split('\t')
	strLG	= cell[0]
	strSC	= cell[1]
	strSTRAND	= cell[2]
	if strSC in done:
		print(strSC)
		continue
	done.append(strSC)
