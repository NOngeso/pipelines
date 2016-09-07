#!/usr/bin/python
from __future__ import print_function
import sys
sys.path.append('/ref/analysis/pipeline/')
import kang

file_cdhit_fa = 'transcripts.fasta.transdecoder.cds.cdhit'
dic_cdhit_fa  = kang.Fasta2dic_all(file_cdhit_fa)

Gene_list       = []
Transcript_list = []

for key in dic_cdhit_fa.keys():
	cell = key.split()
	Gene_list.append(cell[1])
	Transcript_list.append(cell[0])

Gene_list       = set(Gene_list)
Transcript_list = set(Transcript_list)




file_gff = 'transcripts.fasta.transdecoder.gff3'
Outfile = open(file_gff+'.cdhit.gff3','w')
for line in open(file_gff):
	if line.strip() == '':
		continue
	cell    = line.strip().split('\t')
	strT    = cell[2]
	info    = cell[-1]
	dicinfo = dict(zip([x.split('=')[0] for x in info.split(';')],[x.split('=')[1] for x in info.split(';')]))
	if strT == 'gene':
		if dicinfo['ID'] in Gene_list:
			print(line.strip(),file=Outfile)
	elif strT == 'mRNA':
		if dicinfo['ID'] in Transcript_list:
			print(line.strip(),file=Outfile)
	else:
		if dicinfo['Parent'] in Transcript_list:
                        print(line.strip(),file=Outfile)
		

