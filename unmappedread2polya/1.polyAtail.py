#!/usr/bin/python

import pandas as pd
import sys
#1.read mapped in proper pair
#2.read unmapped
#3.mate unmapped
#4.read reverse strand
#5.mate reverse strand
#6.first in pair
#7.second in pair
#8.not primary alignment
#9.read fails platform/vendor quality checks
#10.read is PCR or optical duplicate
#11.supplementary alignment
def flagparser(intFlag):
    dic_key = ['supplementary alignment','read is PCR or optical duplicate','read fails platform/vendor quality checks',\
               'not primary alignment','second in pair','first in pair','mate reverse strand','read reverse strand','mate unmapped',\
              'read unmapped','read mapped in proper pair','read prraired']
    bFlag     = "{0:b}".format(intFlag)
    addzero   = '0'*(12-len(bFlag))
    dic_value = list(addzero+bFlag)
    #print (addzero+bFlag)
    return dict(zip(dic_key,dic_value))

file_polyT = 'polyTfirstreads.txt'
list_polyT = open(file_polyT).readlines()
list_polyT = set([x.strip() for x in list_polyT])

for line in sys.stdin:
    cell    = line.strip().split('\t')
    if cell[0] in list_polyT:
        pass
    else: continue
#    print cell
    intFlag = int(cell[1])
    dic     = flagparser(intFlag)
    seq     = cell[9]
    if dic['read unmapped'] == '1' and dic['second in pair'] == '1':
        #print 1
        if seq[-10:] == 'A'*10:
            print line.strip()
        

