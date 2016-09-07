#!/usr/bin/python

'''
chromosome_2    psl2gff exon    6078249 6079319 .       -       .       TR12395|c0_g4_i1|m.19660
chromosome_2    psl2gff exon    6079511 6079626 .       -       .       TR12395|c0_g4_i1|m.19660
chromosome_2    psl2gff exon    6080032 6080557 .       -       .       TR12395|c0_g4_i1|m.19660
chromosome_6    psl2gff exon    3843396 3843491 .       -       .       TR51629|c0_g1_i1|m.94736
chromosome_6    psl2gff exon    3843757 3843840 .       -       .       TR51629|c0_g1_i1|m.94736
chromosome_6    psl2gff exon    3844092 3844172 .       -       .       TR51629|c0_g1_i1|m.94736
chromosome_6    psl2gff exon    3844466 3844564 .       -       .       TR51629|c0_g1_i1|m.94736
chromosome_6    psl2gff exon    3844837 3844920 .       -       .       TR51629|c0_g1_i1|m.94736
chromosome_6    psl2gff exon    3845164 3845245 .       -       .       TR51629|c0_g1_i1|m.94736
chromosome_6    psl2gff exon    3845344 3845409 .       -       .       TR51629|c0_g1_i1|m.94736
chromosome_6    psl2gff exon    3845519 3845586 .       -       .       TR51629|c0_g1_i1|m.94736
chromosome_12   psl2gff exon    4127815 4130100 .       +       .       TR28522|c0_g1_i3|m.53384
chromosome_12   psl2gff exon    4130307 4130618 .       +       .       TR28522|c0_g1_i3|m.53384
chromosome_1    psl2gff exon    3743395 3744474 .       +       .       TR18418|c1_g2_i1|m.31424
'''
from __future__ import print_function
import pandas as pd
import sys

file_gff  = sys.argv[1]#'Trinity.fasta.transdecoder.mRNA.complete.blat.Creinhardtii_281_v5.0.fa.psl.cp.tophit.gff'
Outfile = open(file_gff+'.mRNA.gff','w')
df_gff    = pd.read_csv(file_gff,sep='\t',header=None)

df_gff_ix = df_gff.set_index(8)
df_gff_ix.sortlevel(inplace=True)
result = []

for i in set(df_gff_ix.index):
	print(i)
        df    = df_gff_ix.loc[i]
        left  = df[3]
        right = df[4]
        try:
                leftmin  = min(left)
                rightmax = max(right)
                chromosome = df[0][0]
                strand     = df[6][0]

        except:
                leftmin = left
                rightmax = right
                chromosome  = df[0]
                strand      = df[6]
        result.append([chromosome, 'psl2gff','mRNA',leftmin,rightmax,'.',strand,'.','ID=%s'%i])
#result = list(set(result))
for i in df_gff.index:
	chromosome = df_gff.loc[i][0]
	left       = df_gff.loc[i][3]
	right      = df_gff.loc[i][4]
	strand     = df_gff.loc[i][6]
	name       = df_gff.loc[i][8]
	#result.append([chromosome, 'psl2gff','exon',left,right,'.',strand,'.','ID=EXON.%s;Parent=%s'%(name,name)])
	result.append([chromosome, 'psl2gff','CDS',left,right,'.',strand,'.','ID=CDS.%s;Parent=%s'%(name,name)])
result.sort(key=lambda x : x[8])
result.sort(key=lambda x : int(x[4]),reverse=True)
result.sort(key=lambda x : int(x[3]))
result.sort(key=lambda x : x[0])
for cell in result :
	print('\t'.join(map(str,cell)),file=Outfile)


