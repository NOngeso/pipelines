#!/usr/bin/python3 

import glob,kang

fa_list = glob.glob('*.prealn.fa')
orthologs = 'orthologs.txt'
dicAHD2seq = kang.Fasta2dic('all.fa')


dicA2B = {}
for line in open(orthologs):
	cell 	= line.strip().split()
	strA 	= cell[0]
	strB	= cell[1]
	try:
		dicA2B[strA].append(strB)
	except KeyError:
		dicA2B[strA] = [strB]
	try:
		dicA2B[strB].append(strA)
	except KeyError:
		dicA2B[strB] = [strA]


for file_fa in fa_list:
	dicHD2seq 	= kang.Fasta2dic(file_fa)
	dicSPCS2GN 	= {}
	for strHD in dicHD2seq:
		spcs 	= strHD.split('|')[0]
		gn	= strHD.split('|')[1]
		dicSPCS2GN[spcs] = strHD
	try:
		add_list 	= dicA2B[dicSPCS2GN['VRA']]
	except KeyError:
		continue
	add_spcs_num 	= len(set([x.split('|')[0] for x in add_list]))
	print(add_list)
	if add_spcs_num == 21:
		pass
	else: continue

	Outfile = open(file_fa+'.good.fa','w')
	for spcs in add_list:
		print('>'+spcs,file=Outfile)
		print(dicAHD2seq[spcs],file=Outfile)
	for strHD in dicHD2seq:
		print('>'+strHD,file=Outfile)
		print(dicHD2seq[strHD],file=Outfile)
	Outfile.close()



		
	
	 
