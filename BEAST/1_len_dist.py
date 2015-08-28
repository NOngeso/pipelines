#!/usr/bin/python3

import glob,kang,numpy,subprocess

file_list = glob.glob('*.prealn.fa')
for file_in in file_list:
	dic = kang.Fasta2dic(file_in)
	len_list = []
	for key in dic:
		len_list.append(len(dic[key]))
	fAvr 	= numpy.average(len_list)
	fStd	= numpy.std(len_list)
	if fStd/fAvr > 0.1:
		continue
	#print(numpy.average(len_list),numpy.std(len_list))
	subprocess.call('cp %s %s.ok.fa'%(file_in,file_in),shell=True)
