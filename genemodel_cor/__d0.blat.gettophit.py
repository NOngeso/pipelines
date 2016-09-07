#!/usr/bin/python
from __future__ import print_function
import sys
file_psl = sys.argv[1] #'transcripts.fasta.transdecoder.cds.psl.cp'

file_stringtie_gff = 'intron3000.merge.sorted.bam.stringtie.gff.strandcor.gff'
'''
chromosome_1    StringTie       transcript      178720  179282  1000    +       .       gene_id "STRG.1"; transcript_id "STRG.1.1"; cov "6.120189"; FPKM "0.054489"; TPM "0.125072";
chromosome_1    StringTie       exon    178720  179282  1000    +       .       gene_id "STRG.1"; transcript_id "STRG.1.1"; exon_number "1"; cov "6.120189";
chromosome_1    StringTie       transcript      179986  180711  1000    -       .       gene_id "STRG.2"; transcript_id "STRG.2.1"; cov "7.936158"; FPKM "0.070657"; TPM "0.162182";
chromosome_1    StringTie       exon    179986  180139  1000    -       .       gene_id "STRG.2"; transcript_id "STRG.2.1"; exon_number "1"; cov "7.402596";
chromosome_1    StringTie       exon    180238  180342  1000    -       .       gene_id "STRG.2"; transcript_id "STRG.2.1"; exon_number "2"; cov "11.794286";
chromosome_1    StringTie       exon    180617  180711  1000    -       .       gene_id "STRG.2"; transcript_id "STRG.2.1"; exon_number "3"; cov "4.536842";
'''
dicStringTie = {}
for line in open(file_stringtie_gff):
	cell = line.strip().split('\t')
	if cell[2] == 'transcript':
		transcriptID = cell[-1].split(';')[1].split()[1].strip('"')
		dicStringTie[transcriptID] = cell[0]

Outfile = open(file_psl+'.tophit','w')
dic = {}
i = 0
for line in open(file_psl):
        cell = line.strip().split('\t')
        try:
		if cell[0] == 'match':
			continue
		name = cell[9]
	except IndexError:
		continue
	if cell[10] != cell[12]:
		continue
	if cell[4] != '0':
		continue
	if cell[5] != '0':
		continue
	putative_chromosome = dicStringTie[name.split('|')[0]]
	if cell[13] != putative_chromosome:
		continue
        try:
                dic[name].append(cell)
        except KeyError:
                dic[name] = [cell]
for name in dic:
        hitlist = dic[name]
	hitlist.sort(key=lambda x : int(x[7]),reverse = False)
        hitlist.sort(key=lambda x : int(x[0]),reverse = True)
        print ('\t'.join(hitlist[0]),file=Outfile)

