#!/usr/bin/python3
# description : calculate N 50 value of assembly in fasta 
# usage       : python3 [thisfile] [fasta]
from __future__ import print_function
import sys,kang

file_fa = sys.argv[1] # fasta

dicHD2Seq = {}
bulk = open(file_fa).read()
bulk_list= bulk.split('>')
for each in bulk_list:
	if each == '':
		continue
	strHD  = each.split('\n')[0]
	strSeq = ''.join(each.split('\n')[1:]).upper().replace('N','')
	if len(strSeq) < 200:
		continue
	dicHD2Seq[strHD] = strSeq

dicHD2Seq_list = list(dicHD2Seq.keys())
dicHD2Seq_list.sort(key = lambda x:len(dicHD2Seq[x]))
total_length = 0
for strHD in dicHD2Seq_list:
	total_length += len(dicHD2Seq[strHD])
print('Total length:',total_length)
print('Number of contig or scaffolds',len(dicHD2Seq_list))
print('Max/min length of contigs or scaffolds',len(dicHD2Seq[dicHD2Seq_list[-1]]),len(dicHD2Seq[dicHD2Seq_list[0]]))
N10 = int(total_length * 0.9)
N20 = int(total_length * 0.8)
N30 = int(total_length * 0.7)
N40 = int(total_length * 0.6)
N50 = int(total_length * 0.5)
N60 = int(total_length * 0.4)
N70 = int(total_length * 0.3)
N80 = int(total_length * 0.2)
N90 = int(total_length * 0.1)
N_list = [N10,N20,N30,N40,N50,N60,N70,N80,N90]
stack = 0
j = 0
print('N##','Length of N##','Number of over N##')
for strHD in dicHD2Seq_list:
	i = 10
	for N in N_list:
		if stack <= N < stack+ len(dicHD2Seq[strHD]):
			print('N%d'%i,len(dicHD2Seq[strHD]),len(dicHD2Seq_list)-j,sep='\t')
		i += 10
	stack += len(dicHD2Seq[strHD])
	j += 1

	

