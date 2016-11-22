#!/usr/bin/python3

import sys

file_non_anchor = sys.argv[1]#'Va.ref.all.gff.maker.gff.sorted.gff'
Prefix		= 'Vradi'
Outfile		= open('%s.scaffold.gff'%Prefix,'w')
Outfile_map	= open('%s.scaffold.map'%Prefix,'w')
dicPGN2NGN = {}
strScaff = ''
for line in open(file_non_anchor):
	cell	= line.strip().split('\t')
	if strScaff != cell[0]:
		i = 10
	strScaff	= cell[0]
	if 'LG' in cell[0]:
		strS	= 'Chr' # state
		intChr = int(cell[0].replace('LG',''))
	elif 'scaffold' in cell[0]:
		strS 	= 'scaffold'
		intChr	= int(cell[0].replace('scaffold_',''))
	elif 'SuperScaf' in cell[0]:
		strS	= 'SuperScaf'
		intChr 	= int(cell[0].replace('SuperScaf_',''))
	strP		= cell[1]
	strT		= cell[2]
	strL1		= cell[3]
	strL2		= cell[4]
	strScore	= cell[5]
	strStrand	= cell[6]
	strFrame	= cell[7]
	strINFO		= cell[8]
	strINFO_list	= cell[8].split(';')
	if strT == 'contig':
		continue
	if strT == 'gene':
		strPGN          = '-'.join(strINFO.split(';')[1].split(':')[0].split('-')).replace('Name=','')
		try:
			strNGN = dicPGN2NGN[strPGN]
		except KeyError:
			if strS == 'scaffold':
				strNGN = '%s%04ds%05d'%(Prefix,intChr,i)
				dicPGN2NGN[strPGN] = strNGN
				i += 10
			elif strS == 'SuperScaf':
				strNGN = '%s%04dS%05d'%(Prefix,intChr,i)
				dicPGN2NGN[strPGN] = strNGN
				i += 10
			elif strS == 'Chr':
				strNGN = '%s%02dg%05d'%(Prefix,intChr,i)
				dicPGN2NGN[strPGN] = strNGN
				i += 10
		print(strScaff,strP,strT,strL1,strL2,strScore,strStrand,strFrame,strINFO_list[0]+';'+'Name='+strNGN+';'+';'.join(strINFO_list[2:]),sep='\t',file=Outfile)
	elif strT == 'mRNA':
		strPGN          = '-'.join(strINFO.split(';')[2].split(':')[0].split('-')[:-2]).replace('Name=','')
		strPTGN		= '-'.join(strINFO.split(';')[2].split(':')[0].split('-')).replace('Name=','')
		strNT   = strINFO.split(';')[2].split(':')[0].split('-')[-1]
#		print(line)
		strNGN  = dicPGN2NGN[strPGN]
		print(strScaff,strP,strT,strL1,strL2,strScore,strStrand,strFrame,strINFO_list[0]+';'+'Name='+strNGN+'.'+strNT+';'+';'.join(strINFO_list[3:]),sep='\t',file=Outfile)
		dicPGN2NGN[strPTGN] = strNGN+'.'+strNT
	elif strT == 'exon':
		#scaffold_0	maker	exon	1530382	1530477	.	+	.	ID=maker-scaffold_0-augustus-gene-15.28-mRNA-1:exon:1063;Parent=maker-scaffold_0-augustus-gene-15.28-mRNA-1,maker-scaffold_0-augustus-gene-15.28-mRNA-2
#		strPGN          = '-'.join(strINFO.split(';')[1].split(':')[0].split('-')[:-2]).replace('Parent=','')
		strPGN_list	= strINFO.split(';')[1].split(',')
		for strPGN_pre in strPGN_list:
			strPGN	= '-'.join(strPGN_pre.replace('Parent=','').split('-')[:-2])
			strNT   = strPGN_pre.replace('Parent=','').split('-')[-1]
	#		print(line)
			strNGN  = dicPGN2NGN[strPGN]
			print(strScaff,strP,strT,strL1,strL2,strScore,strStrand,strFrame,strINFO_list[0]+';'+'Parent='+strNGN+'.'+strNT+';'+';'.join(strINFO_list[2:]),sep='\t',file=Outfile)
	else:
		strPGN          = '-'.join(strINFO.split(';')[1].split(':')[0].split('-')[:-2]).replace('Parent=','')
		strNT	= strINFO.split(';')[1].split(':')[0].split('-')[-1]
	#	print(line)
		strNGN	= dicPGN2NGN[strPGN]
		print(strScaff,strP,strT,strL1,strL2,strScore,strStrand,strFrame,strINFO_list[0]+';'+'Parent='+strNGN+'.'+strNT+';'+';'.join(strINFO_list[2:]),sep='\t',file=Outfile)

for strPGN in dicPGN2NGN:
	print(strPGN,dicPGN2NGN[strPGN],sep='\t',file=Outfile_map)
