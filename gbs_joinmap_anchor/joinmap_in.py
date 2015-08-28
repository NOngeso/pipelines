#!/usr/bin/python3
import sys
file_in = sys.argv[1] #'test.gt'
Outfile = open(file_in+'.jm','w')
for line in open(file_in):
	cell = line.strip().split('\t')
	strLocus	= cell[0].split('_')[0].replace('caffold','').replace('SuperScaf','ss')+cell[0].split('_')[1]+'_'+str(int(cell[0].split('_')[2]) + int(cell[1]))
	strRef		= cell[2]
	strGT_list	= cell[3].split(',')
	strJoinmapGT_list  = []
	for strGT in strGT_list:
		if strGT == strRef:
			strJoinmapGT_list.append('a')
		elif strGT == 'N':
			strJoinmapGT_list.append('-')		
		elif strGT == 'H':
			strJoinmapGT_list.append('h')
		else:
			strJoinmapGT_list.append('b')
	print(strLocus,'\t'.join(strJoinmapGT_list),sep='\t',file=Outfile)
