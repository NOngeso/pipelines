#!/usr/bin/python3

import sys

file_in = sys.argv[1] # *.kaks

spcs = sys.argv[2] # species name ex) VAG

file_gff = sys.argv[3] #'Va2Vr.gff'
Outfile = open(file_in+'.Rin','w')
dicGN2L	= {}
for line in open(file_gff):
	cell = line.strip().split('\t')
	strChr	= cell[0]
	strGN	= cell[1]
	strL	= cell[2]
	dicGN2L[strGN] = strChr, strL  

for line in open(file_in):
	if line[0] == '#':
		continue
	cell = line.strip().split('\t')
	strSB 	= cell[0]
	strGN1	= cell[1]
	strGN2	= cell[2]
	if 'Vang' in strGN1:
		strTGN = strGN1
		strSGN = strGN2
	else: 
		strTGN = strGN2
		strSGN = strGN1
	strKs	= cell[5]
	strTChr, strTL	= dicGN2L[strTGN]
	strSChr, strSL	= dicGN2L[strSGN]
	try:
		print(spcs,strSB, strTChr, strTL,str(int(strSChr.replace('Chr','').replace('Vr','').replace('Vr',''))), strSL, strKs,sep='\t',file=Outfile)	
	except ValueError:
		pass
