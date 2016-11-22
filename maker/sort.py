#!/usr/bin/python3 

import sys

file_in = sys.argv[1]
Outfile = open(file_in+'.sorted.gff3','w')

#########################################################################
list_in = open(file_in).readlines()
list_in.sort(key=lambda x : int(x.split('\t')[3]))
list_in.sort(key=lambda x : x.split('\t')[0])

for line in list_in:
	print(line.strip(),file=Outfile)


