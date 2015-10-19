#!/usr/bin/python3 
import kang
file_fa = 'final.assembly_40k_GSFLX_pseudo5k.fasta.nonATGC.fa'
file_in = 'scaffolds.txt'
Outfile = open('superscaf.fa','w')

dicHD2seq = kang.Fasta2dic(file_fa)
done = []
i = 1
for line in open(file_in):
	cell = line.strip().split('\t')
	strLink = 'SuperScaf_%d'%i
	i += 1
	seq = ''
	for ecell in cell:
		strSC  = ecell.split(',')[0].replace('s','scaffold_') # scaffold name
		if strSC in done:
			print(strSC)
			print('!!!')
			exit()
		done.append(strSC)
		strOrt = ecell.split(',')[1][0] # Orientation F or R
		if strOrt == 'F':
			if seq == '':
				seq += dicHD2seq[strSC]
			else:
				seq += 'N'*500 + dicHD2seq[strSC]
		else : 
			if seq == '':
				seq += dicHD2seq[strSC]
			else:
				seq += 'N'*500 + kang.rev_comp(dicHD2seq[strSC])
	print('>'+strLink,file=Outfile)
	print(seq,file=Outfile)
for strHD in dicHD2seq:
	if strHD in done:
		continue
	print('>'+strHD,file=Outfile)
	print(dicHD2seq[strHD],file=Outfile)


		
		
