#!/usr/bin/python3

import kang,sys

file_joo	= sys.argv[1]#'joinmap.ml.n10.txt_ordered.txt'
dicHD2Seq 	= kang.Fasta2dic('superscaf.fa')
dicLG2Seq 	= {}
LGIncludedSC	= []
Outfile_chr	= open(file_joo+'_'+'Pseudo_chr.fa','w')
Outfile_scaff	= open(file_joo+'_'+'non_anchored_scaffolds.fa','w')
for line in open(file_joo):
	if line[0] == '#' or line.strip() == '':
		continue
	cell 	= line.strip().split('\t')
	strLG 	= cell[0]
	print(cell)
	strSC	= cell[1].replace('*','')
	if 'SS' in strSC:
		strSC = strSC.replace('SS','SuperScaf_')
	else : strSC = strSC.replace('s','scaffold_') 
	LGIncludedSC.append(strSC)
	strD	= cell[2]
	if strD == 'F':
		strSeq = dicHD2Seq[strSC]
	elif strD == 'R':
		strSeq = kang.rev_comp(dicHD2Seq[strSC])
	else : strSeq = dicHD2Seq[strSC]
	try:
		dicLG2Seq[strLG] += 'N'*500+strSeq
	except KeyError:
		dicLG2Seq[strLG] = strSeq
for strLG in dicLG2Seq:
	print('>'+strLG,file=Outfile_chr)
	print(dicLG2Seq[strLG],file=Outfile_chr)
for strHD in dicHD2Seq:
	if strHD in LGIncludedSC:
		continue
	print('>'+strHD,file=Outfile_scaff)
	print(dicHD2Seq[strHD],file=Outfile_scaff)

	
	
	


