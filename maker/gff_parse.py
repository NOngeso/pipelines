#!/usr/bin/python3 
import sys
file_gff = sys.argv[1] #redbean.ref.all.gff.maker.gff
gff_list = open(file_gff).readlines()
Outfile = open(file_gff+'.sorted.gff','w')
gff_list.sort(key = lambda x:int(x.split('\t')[3]))
gff_list.sort(key = lambda x:x.split('\t')[0])
for line in gff_list:
	cell = line.strip().split('\t')
	strChr		= cell[0]
	strMaker	= cell[1]
	strType		= cell[2]
	strL1		= cell[3]
	strL2		= cell[4]
	strE1		= cell[5]
	strStrand	= cell[6]
	strE2	 	= cell[7]
	strINFO_list	= cell[8].split(';')
	strINFO		= strINFO_list
	
	
	print(line.strip(),file=Outfile)
