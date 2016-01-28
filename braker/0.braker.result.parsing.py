#!/usr/bin/python
from __future__ import print_function
import sys
#file_gff = 'alyrata20160127.augustus.gff'
file_gff = sys.argv[1]
Outfile = open(file_gff+'.scored.txt','w')
for line in open(file_gff):
	if line[0] != '#':
		continue
	if line[0:7] == '# start':
		genename = line.strip().replace('# start gene ','')
	if line[0:len('# % of transcript supported by hints (any source):')] == '# % of transcript supported by hints (any source):':
		evidence_value = line.strip().replace('# % of transcript supported by hints (any source): ','')
		print(genename, evidence_value,sep='\t',file=Outfile)

