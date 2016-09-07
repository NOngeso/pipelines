#!/usr/bin/python
from __future__ import print_function


file_gff = 'intron3000.merge.sorted.bam.stringtie.gff'
Outfile = open(file_gff+'.strandcor.gff','w')
for line in open(file_gff):
	if line[0] == '#':
		continue
	cell = line.strip().split('\t')
	strand = cell[6]
	if strand == '+':
		cell[6] = '-'
	else:
		cell[6] = '+'
	print('\t'.join(cell),file=Outfile)
	
