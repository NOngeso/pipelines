#!/usr/bin/python3 

import subprocess,glob,sys

bcf_list = glob.glob('*.bcf')

genotype_list = []

strChr 	= sys.argv[1]
strL	= sys.argv[2]
Outfile = open(strChr+'_'+strL+'.gt.out','w')
for bcf in bcf_list:
#	print('bcftools view %s %s:%s-%s |  grep -v "#"'%(bcf,strChr,strL,strL))
	result 	= subprocess.Popen('bcftools view %s %s:%s-%s |  grep -v "INDEL" |grep -v "#"'%(bcf,strChr,strL,strL),shell=True,stdout=subprocess.PIPE).communicate()[0].decode().split('\t')
	if result[0] == '':
		genotype_list.append('N')
		continue
	ref 	= result[3]
	var	= result[4].split(',')[0]
	hetero = result[9].split(':')[0]
	fQuality = float(result[5])
	strInfo_list = result[7].split(';')
	dicInfo = {}
	for strInfo in strInfo_list:
		try:
			key = strInfo.split('=')[0]
			value = strInfo.split('=')[1]
			dicInfo[key] = value
		except IndexError:
			pass
#	print(dicInfo)
	try:
		iDepth	= int(dicInfo['DP'])
		iMQ	= int(dicInfo['MQ'])
	except KeyError:
		genotype_list.append('N')
		continue
#	print(fQuality,iDepth,iMQ)
	if fQuality > 5 and iDepth > 10 and iMQ > 3:
		pass
	else: 
		genotype_list.append('N')
		continue
	if hetero == '0/1':
		genotype_list.append('H')
		continue
	if var == '.':
		genotype_list.append(ref)
	else: genotype_list.append(var) 

print(strChr,strL,ref,','.join(genotype_list),sep='\t',file=Outfile)
Outfile.close()
