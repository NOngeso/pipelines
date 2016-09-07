#!/usr/bin/python
from __future__ import print_function
file_cds_tophit = 'transcripts.fasta.transdecoder.cds.cdhit.psl.tophit'
file_mRNA = 'transcripts.fasta.transdecoder.mRNA.psl.tophit'
dic = {}
Outfile = open(file_mRNA+'.cdhitonly','w')


for line in open(file_cds_tophit):
	cell = line.strip().split('\t')
	genename = cell[9]
	dic[genename] = 1

for line in open(file_mRNA):
	cell = line.strip().split('\t')
	genename = cell[9]
	try:
		if dic[genename]:
			print(line.strip(),file=Outfile)
	except KeyError:
		pass

