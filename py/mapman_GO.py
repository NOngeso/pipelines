#!/usr/bin/python3
import codecs
file_Gmax = 'Gmax_189.txt'

dicGM2VA = {}
file_match	= 'Gm2Va.collinearity.kaks.recentpeakWGD'
Outfile	= 	open('Va.GO.txt','w')
for line in open(file_match):
	if line[0] == '#':
		continue
	cell 	= line.strip().split('\t')
	if 'Va' in cell[1]:
		strVA 	= cell[1]
		strGM	= cell[2].upper().split('.')[0]
	else:
		strVA   = cell[2]	
		strGM   = cell[1].upper().split('.')[0]
	fKs	= float(cell[3])
	try: 
		if dicGM2VA[strGM]:
			if dicGM2VA[strGM][0] == strVA:
				pass
			elif dicGM2VA[strGM][1] < fKs:
				pass
			else:
#				print(strGM,dicGM2VA[strGM], strVA,fKs)
				dicGM2VA[strGM] = [strVA,fKs]
	except KeyError:
		dicGM2VA[strGM] = [strVA,fKs]

for line in open(file_Gmax,encoding='utf-8'):
	if line[0:7] == 'BINCODE':
		continue
	cell = line.strip().split('\t')
	strGN = cell[2].strip("'")
	if strGN[0:5] == 'glyma':
		pass
	else: continue
	try:
		strVa = dicGM2VA[strGN.upper().split('|')[0].split('.')[0]][0]
		print(cell[0],cell[1],"'"+strVa+"'",'\t'.join(cell[3:]),sep='\t',file=Outfile)
	except KeyError:
		pass
