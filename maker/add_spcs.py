#!/usr/bin/python3 

import sys
Outfile = open(sys.argv[2]+'.Rin','w')
for line in open(sys.argv[1]):
	print(sys.argv[2],line.strip(),sep=',',file=Outfile)
