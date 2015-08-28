#!/usr/bin/python3 

file_in = 'joinmap.rg.n20.txt_ordered.txt'
Outfile = open(file_in+'.sort','w')
list_in_pre = open(file_in).readlines()
list_in = [ x for x in list_in_pre if x[0] != '#'] 
list_in.sort(key=lambda x:float(x.split('\t')[3]))
list_in.sort(key=lambda x: int(x.split('\t')[0].replace('LG','')) )
print(''.join(list_in),file=Outfile)
