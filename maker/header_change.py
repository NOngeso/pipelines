#!/usr/bin/python3 
import kang,sys
file_in = sys.argv[1] # *.fasta
file_map	= sys.argv[2] #Vang.scaffold.map 
Outfile = open(file_in+'.hc.fa','w')
dicO2N	= {} #old 2 new
for line in open(file_map):
	cell = line.strip().split('\t')
	dicO2N[cell[0]] = cell[1]
dicHD2seq = kang.Fasta2dic(file_in)
for strHD in dicHD2seq:
	print('>'+dicO2N[strHD],file=Outfile)
	print(dicHD2seq[strHD],file=Outfile)
