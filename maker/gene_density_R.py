#!/usr/bin/python3 
import sys
file_gff = sys.argv[1]
outfile = open(file_gff+'.gd.Rin','w')
for line in open(file_gff):
	cell = line.strip().split('\t')
	strChr = cell[0] 
	strTP	= cell[2]
	strL	= cell[3]
	if 'caf' not in strChr and strTP == 'gene':
		print(strChr, strL,file=outfile,sep='\t')	
