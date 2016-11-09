#!/usr/bin/python
from __future__ import print_function
'''
GTF format
1. seqname - name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix. Important note: the seqname must be one used within Ensembl, i.e. a standard chromosome name or an Ensembl identifier such as a scaffold ID, without any additional content such as species or assembly. See the example GFF output below.
2. source - name of the program that generated this feature, or the data source (database or project name)
3. feature - feature type name, e.g. Gene, Variation, Similarity
4. start - Start position of the feature, with sequence numbering starting at 1.
5. end - End position of the feature, with sequence numbering starting at 1.
6. score - A floating point value.
7. strand - defined as + (forward) or - (reverse).
8. frame - One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
9. attribute - A semicolon-separated list of tag-value pairs, providing additional information about each feature.
'''
import sys

file_gff = 'Creinhardtii_281_v5.5.gene_exons.gff3'
Outfile  = open(file_gff + '.exon.gtf','w')
dicPACID2mrna_geneid = {}

mRNA            = []
exon            = []
five_prime_UTR  = []
CDS             = []
three_prime_UTR = []

for line in open(file_gff):
    if line[0] == '#':
        continue
    cell = line.strip().split('\t')
    strT = cell[2]
    if strT == 'mRNA':
        mRNA.append(line.strip())
    elif strT == 'CDS':
        CDS.append(line.strip())
    elif strT == 'exon':
        exon.append(line.strip())
    elif strT == 'five_prime_UTR':
        five_prime_UTR.append(line.strip())
    elif strT == 'three_prime_UTR':
        three_prime_UTR.append(line.strip())

for line in mRNA:
    cell    = line.split('\t')
    info    = cell[-1]
    dicinfo = dict(zip([x.split('=')[0] for x in info.split(';')],[x.split('=')[1] for x in info.split(';')]))
    Name    = dicinfo['Name']
    PACID   = dicinfo['pacid']
    Parent  = dicinfo['Parent']
    dicPACID2mrna_geneid[PACID] = Name,Parent

for line in exon:
    cell   = line.split('\t')
    strC   = cell[0] # chromosome
    strS   = cell[1] # source
    strF   = cell[2] # Feature
    strSt  = cell[3] # start
    strEd  = cell[4] # end
    strSc  = cell[5] # score
    strSd  = cell[6] # strand
    strFm  = cell[7] # frame
    info   = cell[8] # Info
    dicinfo = dict(zip([x.split('=')[0] for x in info.split(';')],[x.split('=')[1] for x in info.split(';')]))
    PACID   = dicinfo['pacid']
    Name,Parent =  dicPACID2mrna_geneid[PACID]
    print(strC,strS,strF,strSt,strEd,strSc,strSd,strFm,'gene_id "%s";'%(Parent),'transcript_id "%s";'%(Name),sep='\t',file=Outfile)

    

    

