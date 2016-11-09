#!/usr/bin/python3 

file_interpro = 'adzuki.ref.fasta.pep.fa.noasterix.fa.tsv'
dicGN2GO = {}
Outfile = open(file_interpro+'.bingo','w')
for line in open(file_interpro):
	cell = line.strip().split('\t')
	strGN = cell[0]
	try:
		strGO = cell[13].split('|')
		if strGO == []:
			continue
		try:
			dicGN2GO[strGN] += strGO
		except KeyError:
			dicGN2GO[strGN] = strGO
	except IndexError:
		pass
for strGN in dicGN2GO:
	for eachGO in set(dicGN2GO[strGN]):
		if eachGO == '':
			continue
		print(strGN,'=',eachGO.replace('GO:',''),file=Outfile)
	
