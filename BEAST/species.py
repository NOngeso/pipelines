#!/usr/bin/python3 

file_in = 'prespecies'

print('taxon','species',sep='\t')
for line in open(file_in):
	print(line.strip().replace('>',''), line.strip().replace('>','').split('|')[0],sep='\t')
