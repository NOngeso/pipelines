#!/usr/bin/python3 

import glob,kang

mcscanout_list = ['Vr2Cc.collinearity.kaks.recentpeakWGD','Vr2Gm.collinearity.kaks.recentpeakWGD','Vr2Wd.collinearity.kaks.recentpeakWGD','Vra2Vrr.collinearity.kaks.recentpeakWGD']

file_GMAA	= 'GMAA.fasta'
file_VRGA	= 'VRPA.fasta'
file_matchGM	= 'match.GM.txt'
file_matchVRG	= 'match.Vrg.txt'
#file_group	= 'groups.primary.txt'
file_fa		= 'all.fas'
Outfile_seed	= open('seed_orthologs.txt','w')


dicHD2seq	= kang.Fasta2dic(file_fa)

dicGMAA2seq 	= kang.Fasta2dic(file_GMAA)
dicVRGA2seq	= kang.Fasta2dic(file_VRGA)

dicGMA2B	= {}
for line in open(file_matchGM):
	cell 	= line.strip().split('\t')
	strA	= cell[0]
	strB	= cell[1]
	dicGMA2B[strA] = strB
dicVRGA2B	= {}
for line in open(file_matchVRG):
	cell	= line.strip().split('\t')
	strA	= cell[0]
	strB	= cell[1]
	dicVRGA2B[strA] = strB

GMAA_list 	= list(dicGMAA2seq.keys())
VRGA_list	= list(dicVRGA2seq.keys())

dicVR2ortho 	= {}
for mcscanout in mcscanout_list:
	for line in open('../'+mcscanout):
		if line[0] == '#':
			continue
		cell 	= line.strip().split('\t')
		strSB	= cell[0]	
		strG1	= cell[1]
		strG2	= cell[2]
		if 'Vradi' in strG1:
			strVR	= strG1
			strT	= strG2
		else: 
			strVR	= strG2
			strT	= strG1
		try:
			dicVR2ortho[strVR].append(strT)
		except KeyError:
			dicVR2ortho[strVR] = [strT]


seed_cluster = {}
for strVR in dicVR2ortho:
	bPass = 0
	for each in dicVR2ortho[strVR]:
		if each in GMAA_list:
			GMAA	= each
			GMAB	= dicGMA2B[each]
			if GMAB in dicVR2ortho[strVR]:
				bPass += 1
		elif each in VRGA_list:
			VRGA	= each
			VRGB	= dicVRGA2B[each]
			if VRGB in dicVR2ortho[strVR]:
				bPass += 1
	if bPass == 2 and len(dicVR2ortho[strVR]) == 6:
		print(strVR,  dicVR2ortho[strVR],[GMAA,VRGA,GMAB,VRGB],file=Outfile_seed)
		seed_cluster[strVR] = ['VRA|'+strVR,'CCA|'+dicVR2ortho[strVR][0],'VRS|'+dicVR2ortho[strVR][3],'GMAA|'+GMAA,'VRGA|'+VRGA,'GMAB|'+GMAB,'VRGB|'+VRGB]
	elif bPass > 2:
		print('!!!!')
	else:
		continue

dicGN2clst = {}
i = 0

def locus_write(seed_list_all,filename):
	Outfile = open(filename+'.prealn.fa','w')
	for gn in seed_list_all:
		print('>'+gn,file=Outfile)
		try:
			print(dicHD2seq[gn.split('|')[1]],file=Outfile)	
		except KeyError:
			print(dicHD2seq[gn],file=Outfile)
	Outfile.close()

i = 0
for strVR in seed_cluster:
	locus_write(seed_cluster[strVR],'locus%i'%i)	
	i += 1

'''

for line in open(file_group):
	cell	= line.strip().split()
	strCLST	= cell[0]
	strGN_list	= cell[1:]
	spcs_list	= [x.split('|')[0] for x in strGN_list]
	strGN_in_list	= [x.split('|')[1] for x in strGN_list]
	if len(set(spcs_list)) == 28 and len(strGN_list) < 36:
		pass
	else: continue
	for strGN in strGN_list:
		spcs 		= strGN.split('|')[0]
		strGN_in 	= strGN.split('|')[1] 
		if spcs == 'VRA':	
			bPass = 1	
			try:
				seed_list_all = ['VRA|'+strGN_in] + ['CCA|'+seed_cluster[strGN_in][0]] + ['VRS|'+seed_cluster[strGN_in][3]] + ['GMAA|'+seed_cluster[strGN_in][-4]] + ['VRGA|'+seed_cluster[strGN_in][-3]] + ['GMAB|'+seed_cluster[strGN_in][-2]] + ['VRGB|'+seed_cluster[strGN_in][-1]]
				seed_list = [seed_cluster[strGN_in][0],seed_cluster[strGN_in][-4].split('.')[0]+'.1']
				print(seed_list_all,file=Outfile_seed)
			except KeyError:
				bPass = 0
				continue
			for each in seed_list:
				if each not in strGN_in_list:
					bPass = 0 
			break
	if bPass == 1:
		locus_write(seed_list_all,strGN_list,'locus%d'%i)
		i += 1			
'''

				



