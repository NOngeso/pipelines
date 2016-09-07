#!/usr/bin/python
#SRR2132433.4524365.1    355     chromosome_1    1       1       100M    =       151     250     GGGAACCAGCTACTAGATGGTTCGATTAGTCTTTCGCCCCTATACCCAAGTCTGAAAAGCGATTTGCACGTCAGCACATCTACGAGCCTCCACCAGAGTT    CCCFFFFFHHHHHJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJHIIJJJJJJJJJJJJJJHHHHFFFFFDDDECDDDDDDDDDDDDDDDDDDDDBD>>   AS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:100        YT:Z:UU NH:i:4  CC:Z:chromosome_14      CP:i:4148565    HI:i:0
from __future__ import print_function
import re
import sys
import numpy as np
import pickle as pk


def parse_cigar(scigar):
	match = re.findall(r'(\d+)(\w)', scigar)
	M = sum([int(x[0]) for x in match if x[1] == 'M'])
	I = sum([int(x[0]) for x in match if x[1] == 'I'])
	N = sum([int(x[0]) for x in match if x[1] == 'N'])
	return M+I+N

def fasta2dic(file_fa):
	dic  = {}
	bulk = open(file_fa).read()
	bulk_list = bulk.split('>')
	for each in bulk_list:
		if each.strip() == '':
			continue
		contig_name 	= each.split('\n')[0]
		seq 		= ''.join(each.split('\n')[1:])
		dic[contig_name] = seq
	return dic	

file_fa     = '/ref/analysis/References/Creinhardtii/Creinhardtii_281_v5.0.fa'
Outfilename = sys.argv[1]


dicHD2seq = fasta2dic(file_fa) 						###<<< dictionary contig 2 sequence
dicHD2array = {}
for strHD in dicHD2seq.keys():
	seq = dicHD2seq[strHD]
	dicHD2array[strHD] = np.zeros(len(seq),dtype=np.int)				###<<< dicHD2array : contig 2 count array



for line in sys.stdin:
	cell 		= line.strip().split('\t')
	iBridge_len 	= int(cell[8])
	#if 300 < abs(iBridge_len) < max_insert_len:
	#	pass
	#else: continue
	strChr 		= cell[2]
	if strChr == '*':
		continue
	intStart 	= int(cell[3])
	strCigar 	= cell[5]
	fraglen = parse_cigar(strCigar)
	if 1:
		intEnd 	= intStart + fraglen -1
		try:
			dicHD2array[strChr][intStart:intEnd] += 1
		except :
			print(line)
			exit()

#mask = dicHD2array['chromosome_1']>0
#print(dicHD2array['chromosome_1'][mask])

pk.dump( dicHD2array, open( "%s.p"%(Outfilename), "wb" ) )
