#!/usr/bin/python3
import sys
file_in = sys.argv[1] # 'all.gt.seg.txt'
Outfile = open(file_in+'.jm','w')
for line in open(file_in):
	cell = line.strip().split('\t')
	strLocus	= cell[0].replace("caffold",'')+'_'+cell[1]
	strRef		= cell[2]
	strGT_list	= cell[3].split(',')
	strJoinmapGT_list  = []
	for strGT in strGT_list:
		if strGT == strRef:
			strJoinmapGT_list.append('a')
		elif strGT == 'N':
			strJoinmapGT_list.append('n')		
		elif strGT == 'H':
			strJoinmapGT_list.append('h')
		else:
			strJoinmapGT_list.append('b')
	print(strLocus,'\t'.join(strJoinmapGT_list),sep='\t',file=Outfile)
