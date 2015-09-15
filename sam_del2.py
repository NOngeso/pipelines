#!/usr/bin/python3

import sys,numpy
import subprocess

file_bam     = sys.argv[1]
Outfile_sam  = open(file_bam.split('/')[-1]+'.sam','w')
Outfile_del  = open(file_bam.split('/')[-1]+'.del','w')
Samtools_out = subprocess.Popen('samtools view %s'%file_bam,shell=True,stdout=subprocess.PIPE).stdout

stack = ''
list_loc1 = []
list_loc2 = []
list_insert_size = []
list_chr = []

	

i = 0
done_chr = []
for line in Samtools_out:
	line    = line.decode('ascii')
	cell    = line.split('\t')
	strChr  = cell[2]
	if strChr not in done_chr:
		print(strChr,"processing")
		done_chr.append(strChr)
	loc1    = int(cell[3])
	loc2    = int(cell[7])
	insert_size = int(cell[8])
	
	if cell[6] != '=': #diff chr mapped
		continue
	if insert_size > 100000: #weird insert size
		continue
	if insert_size < 1000: #normal insert size
		continue 
	if list_loc1 == []:
		list_chr.append(strChr)
		list_loc1.append(loc1)
		list_loc2.append(loc2)
		list_insert_size.append(insert_size)
		continue
	if list_chr[-1] != strChr:
		fltIS_avr = numpy.average(list_insert_size)
		fltIS_std = numpy.std(list_insert_size)
		if len(list_loc1) > 10 and fltIS_avr * 0.1 > fltIS_std:
			print(strChr,min(list_loc1),max(list_loc1),min(list_loc2),max(list_loc2),fltIS_avr,fltIS_std,len(list_loc1),file=Outfile_del,sep='\t')
			print(stack.strip(),file=Outfile_sam)
		stack = ''
		list_loc1 = []
		list_loc2 = []
		list_insert_size = []
		continue
	
	if abs(list_loc1[-1] - loc1) > 10000: #reads gap supporting same del : 1000
		#print line
		#print abs(list_loc1[-1] - loc1)
		fltIS_avr = numpy.average(list_insert_size)
		fltIS_std = numpy.std(list_insert_size)
		if len(list_loc1) > 10 and fltIS_avr * 0.1 > fltIS_std:
			print(strChr,min(list_loc1),max(list_loc1),min(list_loc2),max(list_loc2),fltIS_avr,fltIS_std,len(list_loc1),file=Outfile_del,sep='\t')
			print(stack.strip(),file=Outfile_sam)
		stack = ''
		list_loc1 = []
		list_loc2 = []
		list_insert_size = []
		continue
	if abs(list_loc2[-1] - loc2) > 10000: #reads gap supporting same del : 1000
		#print line
		#print abs(list_loc2[-1] - loc2)
		fltIS_avr = numpy.average(list_insert_size)
		fltIS_std = numpy.std(list_insert_size)
		if len(list_loc1) > 10 and fltIS_avr * 0.1 > fltIS_std:
			print(strChr,min(list_loc1),max(list_loc1),min(list_loc2),max(list_loc2),fltIS_avr,fltIS_std,len(list_loc1),file=Outfile_del,sep='\t')
			print(stack.strip(),file=Outfile_sam)
		stack = ''
		list_loc1 = []
		list_loc2 = []
		list_insert_size = []
		continue
	
	list_chr.append(strChr)
	list_loc1.append(loc1)
	list_loc2.append(loc2)
	list_insert_size.append(insert_size)
	stack += line

		
