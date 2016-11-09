#!/usr/bin/python3 

import subprocess,glob,kang,os

file_orth 	= 'orthologs.txt'
Outfile 	= open(file_orth+'.shared','w')
dicHD2seq 	= kang.Fasta2dic('../cds/all.cds.fa')
dicA2B 		= {}
for line in open(file_orth):
	cell 	= line.strip().split('\t')
	strA 	= cell[0]
	strB	= cell[1]
	if strA.split('|')[0] == 'VAG':
		try:
			dicA2B[strA].append(strB)
		except KeyError:
			dicA2B[strA] = [strA,strB]
	else: 
		try:
			dicA2B[strB].append(strA)
		except KeyError:
			dicA2B[strB] = [strB,strA]


def convert_axt(in_file,out_file):
	Outfile = open(out_file,'w')
	dicHD2seq = kang.Fasta2dic(in_file)
	print(out_file,file=Outfile)
	for strHD in dicHD2seq:
		print(dicHD2seq[strHD],file=Outfile)
	Outfile.close()

def get_div(fKs):
	constant        = 5.17 * 10**(-3)
	divergence_time = fKs/(2*constant)
	return(divergence_time)

def get_Ks_list(listin):
	i 	= 1
	dic	= {}
	done 	= []
	listin.sort(key=lambda x : x.split('|')[0])
	for G1 in listin:
		for G2 in listin[i:]:
			G1_spcs = G1.split('|')[0]
			G2_spcs = G2.split('|')[0]
			Outfile = open('./temp/%s.%s.temp.fa'%(G1_spcs,G2_spcs),'w')
			if [G1_spcs,G2_spcs] in done:
				continue
			else:
				done.append([G1_spcs,G2_spcs])
			print('>'+G1_spcs,file=Outfile)
			try:
				print(dicHD2seq[G1],file=Outfile)
			except KeyError:
				print(dicHD2seq[G1.replace('_1','')],file=Outfile)
			print('>'+G2_spcs,file=Outfile)
			try:
				print(dicHD2seq[G2],file=Outfile)
			except KeyError:
				print(dicHD2seq[G2.replace('_1','')],file=Outfile)
			Outfile.close()
		i += 1
	subprocess.call('parallel -j 10 /data2/k821209/programs/prank/bin/prank -translate -d={1} -o={1}.aln ::: ./temp/*.temp.fa',shell=True)
	aln_list = glob.glob('./temp/*.best.nuc.fas')
	for alnfile in aln_list:
		convert_axt(alnfile,alnfile+'.axt')
	axt_list = glob.glob('./temp/*.axt')
	#for axtfile in axt_list:
	subprocess.call('parallel -j 10 /data2/k821209/programs/KaKs_Calculator1.2/src/KaKs_Calculator -i {1} -o {1}.kaksout -m LWL > /dev/null ::: ./temp/*.axt',shell=True)
	kaksout_list = glob.glob('./temp/*.kaksout')
	for kaksout in kaksout_list:
		result = open(kaksout).read()
		strKey	= ','.join(result.split('\n')[1].split('\t')[0].split('/')[-1].split('.')[0:2])
		strKs   = result.split('\n')[1].split('\t')[3]
		strKa   = result.split('\n')[1].split('\t')[2]
		strKaKs = result.split('\n')[1].split('\t')[4]
		strP    = result.split('\n')[1].split('\t')[5]
		if strKs == 'NA':
			continue
		if strKa == 'NA':
			continue
		try:
			if dic[strKey]:
				print('!!!!')
				print(strKey)
				print(dic)
				print(listin)
				exit()
		except KeyError:
			dic[strKey] = float(strKs),float(strKa)
	os.system('rm ./temp/*.temp.*')
	return(dic)

def get_peak(Ks_list):
	maxKs           = 2
	fWindow         = 0.01
	dicW2C          = {} # Window 2 count
	for i in range(int(maxKs/fWindow)):
		dicW2C[(fWindow*i,fWindow*(i+1))] = 0
	for strKs in Ks_list:
		for fL, fR in dicW2C:
			if fL < float(strKs) <= fR:
				dicW2C[(fL,fR)] += 1
	dicW2C_list = list(dicW2C.keys())
	dicW2C_list.sort(key = lambda x: dicW2C[x], reverse=True)
	return(dicW2C_list[0])


dicP2Ks = {}
dicP2Ka	= {}
i = 0
for strA in dicA2B:
	listin = dicA2B[strA]
	spcs_pre = [x.split('|')[0] for x in listin]	
	spcs = set(spcs_pre)
	if len(spcs) == 2:
		print(1)
		print('\t'.join(listin),file=Outfile)
		done = []
		listin_filtered = []
		for each in listin:
			eachspcs = each.split('|')[0]
			if eachspcs in done:
				continue
			else:
				done.append(eachspcs)
				listin_filtered.append(each)
		dic = get_Ks_list(listin_filtered)
		#print(dic)
		for key in dic:
			try:
				dicP2Ks[key].append(dic[key][0])
				dicP2Ka[key].append(dic[key][1])
			except KeyError:
				dicP2Ks[key] = [dic[key][0]]
				dicP2Ka[key] = [dic[key][1]]
		print('processing',strA)
#		if i == 10:
#			break
		i += 1

Outfile_ks = open('Ks_dist.txt','w')
for pair in dicP2Ks:
	Ks_list = dicP2Ks[pair]
	peak = get_peak(Ks_list)
	print(pair,get_div((peak[0]+peak[1])/2),','.join(map(str,peak)),','.join(map(str,Ks_list)),sep='\t',file=Outfile_ks)
Outfile_ks.close()
Outfile_ka = open('Ka_dist.txt','w')
for pair in dicP2Ka:
	Ka_list = dicP2Ka[pair]
	peak = get_peak(Ka_list)
	print(pair,get_div((peak[0]+peak[1])/2),','.join(map(str,peak)),','.join(map(str,Ka_list)),sep='\t',file=Outfile_ka)
Outfile_ka.close()

	
