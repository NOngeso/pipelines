#!/usr/bin/python
from __future__ import print_function
file_gff = 'transcripts.fasta.transdecoder.genome.cdhit.gff3'

GeneID = []
Outfile = open(file_gff+'.genecor.gff3','w')
for line in open(file_gff):
	if line[0] == '#' or line.strip() == '':
		continue
	cell = line.strip().split('\t')
	strT = cell[2] 
	if strT == 'gene':
		eGeneID = cell[-1].split(';')[0].replace('ID=','')
		if eGeneID in set(GeneID):
			pass
		else:
			print(line.strip(),file=Outfile)
			GeneID.append(eGeneID)	
	else:
		print(line.strip(),file=Outfile)


