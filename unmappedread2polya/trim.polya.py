#!/usr/bin/python
from __future__ import print_function
file_sam = 'polyareads.txt'
Outfile  = open(file_sam+'.fq','w')
for line in open(file_sam):
    cell = line.strip().split('\t')
    readid = cell[0]
    seq  = cell[9].rstrip('A')
    if len(seq) < 50:
        continue
    qual = cell[10][0:len(seq)]
    print('@'+readid,file=Outfile)
    print(seq,file=Outfile)
    print('+',file=Outfile)
    print(qual,file=Outfile)

      
