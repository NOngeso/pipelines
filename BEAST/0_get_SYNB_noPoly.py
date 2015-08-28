#!/usr/bin/python3 

import glob,kang,glob

#mcscanout_list = ['Vr2Cc.collinearity.kaks.recentpeakWGD','Vr2Gm.collinearity.kaks.recentpeakWGD','Vr2Wd.collinearity.kaks.recentpeakWGD','Vra2Vrr.collinearity.kaks.recentpeakWGD']
#mcscanout_list_pre = 'At2Al.collinearity.kaks.recentpeakWGD.add  At2Br.collinearity.kaks.recentpeakWGD.add  At2Th.collinearity.kaks.recentpeakWGD.add'
#mcscanout_list = mcscanout_list_pre.split()
mcscanout_list = glob.glob('*.add')

file_fa		= 'all.fas'
Outfile_seed	= open('seed_orthologs.txt','w')


dicHD2seq	= kang.Fasta2dic(file_fa)

dicNest2ortho 	= {}
for mcscanout in mcscanout_list:
	for line in open('./'+mcscanout):
		if line[0] == '#':
			continue
		cell 	= line.strip().split('\t')
		strSB	= cell[0]	
		strG1	= cell[1]
		strG2	= cell[2]
		if 'VAG' == strG1[0:3]:
			strNest	= strG1
			strT	= strG2
		else: 
			strNest	= strG2
			strT	= strG1
		try:
			dicNest2ortho[strNest].append(strT)
		except KeyError:
			dicNest2ortho[strNest] = [strT]


seed_cluster = {}
for strNest in dicNest2ortho:
	spcs_list = [x.split('|')[0] for x in dicNest2ortho[strNest]]
	spcs_set = list(set(spcs_list))
	print(spcs_set)
	if len(spcs_set) == len(mcscanout_list) == len(spcs_list):
		print(strNest,  dicNest2ortho[strNest],file=Outfile_seed)
		#seed_cluster[strNest] = ['ATH|'+strNest,'CCA|'+dicNest2ortho[strNest][0],'VRS|'+dicNest2ortho[strNest][3]]
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
			print(dicHD2seq[gn.split('|')[1]],file=Outfile)
	Outfile.close()

i = 0
for strNest in dicNest2ortho:
	spcs_list = [x.split('|')[0] for x in dicNest2ortho[strNest]]
	spcs_set = list(set(spcs_list))
	if len(spcs_set) == len(mcscanout_list) == len(spcs_list):
		locus_write(dicNest2ortho[strNest]+[strNest],'locus%i'%i)	
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

				



