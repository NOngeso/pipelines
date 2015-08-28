#!/usr/bin/python3
#Gm01	39	.	T	.	90	.	DP=26;AF1=0;AC1=0;DP4=18,2,0,0;MQ=36;FQ=-87	PL	0
#Gm01	30	.	T	G	23	.	DP=14;VDB=0.0008;AF1=0.5;AC1=1;DP4=3,0,6,0;MQ=41;FQ=26;PV4=1,0.026,4e-05,1	GT:PL:GQ	0/1:53,0,71:56
#Gm01	79	.	GTT	.	9.1	.	INDEL;DP=44;VDB=0.0124;AF1=0.4997;AC1=1;DP4=16,8,3,2;MQ=39;FQ=-6.6;PV4=1,1,0.08,0.00066	PL	29
import sys
from datetime import datetime

Outfile = open(sys.argv[1]+'.out.fa','w')
start_base = 1
seq = ''
strChr = ''
minD 	= 2
minMQ	= 30
maxD 	= 50
for line in sys.stdin:
	if strChr == '':
		print('done',str(datetime.now()),strChr)
	if line[0] == '#':
		continue
	cell 	= line.strip().split('\t')
	#########################################################################
	if strChr != '' and strChr != cell[0]:
		print('done',str(datetime.now()),strChr)
		print('>'+strChr,file=Outfile)
		print(seq,file=Outfile)
		seq = ''
		start_base = 1
	#########################################################################
	if 'INDEL' == cell[7][0:5]:
		continue # Ignore INDEL sites
	strChr 	= cell[0]
	intLoc	= int(cell[1])
	strRef	= cell[3]
	strVar	= cell[4].split(',')[0]
	dicValues = {}
#	print(line)
	for each in cell[7].split(';'):
		if each == '':
			continue
		key = each.split('=')[0]
		value = each.split('=')[1]
		dicValues[key] = value
	try:	
		intD    = int(dicValues['DP'])
		intMQ	= int(dicValues['MQ'])
	except KeyError:
		continue
	if minD <= intD <= maxD:
		pass
	else: continue

	if intMQ < minMQ:
		continue
	#########################################################################
	while intLoc > start_base:
		seq += 'n'
		start_base += 1
	if intLoc == start_base:
		if strVar == '.':
			seq += strRef
			start_base += 1
		else : 
			seq += strVar
			start_base += 1

	else: 
		print(intLoc, start_base)
		exit()
print('done',str(datetime.now()),strChr)
print('>'+strChr,file=Outfile)
print(seq,file=Outfile)

	
	
	
