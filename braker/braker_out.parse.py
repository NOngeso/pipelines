#!/usr/bin/python3


'''
chromosome_1    AUGUSTUS        gene    23183   29577   0.94    +       .       g1
chromosome_1    AUGUSTUS        transcript      23183   29577   0.94    +       .       g1.t1
chromosome_1    AUGUSTUS        start_codon     23183   23185   .       +       0       transcript_id "g1.t1"; gene_id "g1";
chromosome_1    AUGUSTUS        intron  23251   23370   0.99    +       .       transcript_id "g1.t1"; gene_id "g1";
chromosome_1    AUGUSTUS        intron  23403   24372   0.94    +       .       transcript_id "g1.t1"; gene_id "g1";
chromosome_1    AUGUSTUS        intron  28106   28290   1       +       .       transcript_id "g1.t1"; gene_id "g1";
chromosome_1    AUGUSTUS        intron  28645   28841   1       +       .       transcript_id "g1.t1"; gene_id "g1";
chromosome_1    AUGUSTUS        intron  29092   29346   1       +       .       transcript_id "g1.t1"; gene_id "g1";
chromosome_1    AUGUSTUS        CDS     23183   23250   0.99    +       0       transcript_id "g1.t1"; gene_id "g1";
chromosome_1    AUGUSTUS        CDS     23371   23402   0.99    +       1       transcript_id "g1.t1"; gene_id "g1";
chromosome_1    AUGUSTUS        CDS     24373   28105   0.94    +       2       transcript_id "g1.t1"; gene_id "g1";
chromosome_1    AUGUSTUS        CDS     28291   28644   1       +       1       transcript_id "g1.t1"; gene_id "g1";
chromosome_1    AUGUSTUS        CDS     28842   29091   1       +       1       transcript_id "g1.t1"; gene_id "g1";
chromosome_1    AUGUSTUS        CDS     29347   29577   1       +       0       transcript_id "g1.t1"; gene_id "g1";
chromosome_1    AUGUSTUS        stop_codon      29575   29577   .       +       0       transcript_id "g1.t1"; gene_id "g1";
'''

file_gff = 'braker.try4.augustus.gff.rmsharp.gff'
Outfile = open(file_gff+'.parsed.gff','w')
for line in open(file_gff):
	cell   = line.strip().split('\t')
	strT   = cell[2]
	if strT == 'gene':
		cell[-1] = 'ID='+cell[-1].strip()
	elif strT == 'transcript':
		cell[2]  = 'mRNA'
		cell[-1] = 'ID='+cell[-1]+';'+'Parent='+cell[-1].split('.')[0]
	elif strT == "5'-UTR":
		cell[2]  = 'five_prime_UTR'
		cell[-1] = 'Parent='+cell[-1].split(';')[0].split()[1].strip('"')
	elif strT == "3'-UTR":
		cell[2]  = 'three_prime_UTR'
		cell[-1] = 'Parent='+cell[-1].split(';')[0].split()[1].strip('"')
	else:
		cell[-1] = 'Parent='+cell[-1].split(';')[0].split()[1].strip('"')
	print('\t'.join(cell),file=Outfile)
