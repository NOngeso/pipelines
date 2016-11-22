#!/usr/bin/python3
import sys
file_ipr	= sys.argv[1] # ipr tsv
file_map        = sys.argv[2] #Vang.scaffold.map 
Outfile = open(file_ipr+'.hc.tsv','w')
dicO2N  = {} #old 2 new
for line in open(file_map):
        cell = line.strip().split('\t')
        dicO2N[cell[0]] = cell[1]
for line in open(file_ipr):
	cell = line.strip().split('\t')
	strGN_list = cell[0].split('|') 
	if len(strGN_list) > 1:
		strGN_changed_list = [dicO2N[x] for x in strGN_list]
		print('|'.join(strGN_changed_list),'\t'.join(cell[1:]),sep='\t',file=Outfile)
	else:
		print(dicO2N[strGN_list[0]],'\t'.join(cell[1:]),sep='\t',file=Outfile)
