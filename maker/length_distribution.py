#!/usr/bin/python3 

import sys

file_gff = sys.argv[1] #*.gff

dicmRNA2CDS = {}
dicT2L	= {}
for line in open(file_gff):
	if line[0] == '#' or line.strip() == '':
		continue
	cell = line.strip().split('\t')
	strChr = cell[0]
	strTP	= cell[2]
	strL1	= cell[3]
	strL2	= cell[4]
	strInfo_list = cell[-1].split(';')
	dicInfo = {}
	for strInfo in strInfo_list:
		try:
			key = strInfo.split('=')[0]
			value = strInfo.split('=')[1]
			dicInfo[key] = value
		except IndexError:
			pass
	if strTP == 'CDS':
		try:
			dicmRNA2CDS[dicInfo['Parent']].append([int(strL1),int(strL2)])
		except KeyError:
			dicmRNA2CDS[dicInfo['Parent']] = [[int(strL1),int(strL2)]]
	try:
		dicT2L[strTP].append(int(strL2)-int(strL1))
	except KeyError:
		dicT2L[strTP] = [int(strL2)-int(strL1)]

for T in dicT2L:
	Outfile = open(file_gff + '_' + T + '.ld.out','w')
	print(T,file=Outfile)
	print('\n'.join(map(str,dicT2L[T])),file=Outfile)
	Outfile.close()

Outfile = open(file_gff + '_' + 'intron' + '.ld.out','w')
print('intron',file=Outfile)
for mRNA in dicmRNA2CDS:
	list_in = dicmRNA2CDS[mRNA]
	list_in.sort(key=lambda x : x[0])
	for i in range(len(list_in)):
		try:
			intron_len = list_in[i+1][0] - list_in[i][1]
		except IndexError:
			continue
		if intron_len < 0 :
			print('!')
			print(mRNA)
			continue
		else: 
			print(intron_len,file=Outfile)

		

		


