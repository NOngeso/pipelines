#!/usr/bin/python3

import sys

file_del   = sys.argv[1]
Outfile    = open(file_del+'.out','w')
file_GFF   = '/data/ref/Gmax_189_gene.gff3'
file_annot = '/data/ref/Gmax_189_annotation_info.txt'
#########################################################################
dicGN2Annot = {}
for line in open(file_annot):
	cell     = line.strip().split('\t')
	strGN    = cell[1].split('.')[0]
	strAnnot = cell[-1]
	try:
		if dicGN2Annot[strGN]:
			pass
	except KeyError:
		dicGN2Annot[strGN] = strAnnot  
#########################################################################
dicChr2info_list = {}
for line in open(file_GFF):
	if line[0] == '#':
		continue
	cell   = line.strip().split('\t')
	strChr  = cell[0]
	strTY   = cell[2]
	strLOC1 = cell[3] 
	strLOC2 = cell[4]
	if strTY == 'gene':
		strGN   = cell[-1].split(';')[0].replace('ID=','')
		try: 
			dicChr2info_list[strChr].append([strLOC1,strLOC2,strGN])
		except KeyError:	
			dicChr2info_list[strChr] = [[strLOC1,strLOC2,strGN]]

#########################################################################

def get_annot(strChr,strLoc1,strLoc2):
	strInfo_list = dicChr2info_list[strChr]
	for strInfo in strInfo_list:
		strLoc1_T = strInfo[0]
		strLoc2_T = strInfo[1]
		strGN     = strInfo[2]
		if int(strLoc1) <= int(strLoc1_T) <= int(strLoc2) or int(strLoc1) <= int(strLoc2_T) <= int(strLoc2):
			return([strChr,strLoc1,strLoc2,strGN,dicGN2Annot[strGN]])
	return([strChr,strLoc1,strLoc2,'non'])
def main():
	for line in open(file_del):
		cell    = line.strip().split('\t')
		strChr  = cell[0]
		if 'Gm' not in strChr:
			continue
		strLoc1 = cell[2]
		strLoc2 = cell[3]
		if int(strLoc2) - int(strLoc1) < 0:
			continue
		print(line.strip(),'\t'.join(get_annot(strChr,strLoc1,strLoc2)),sep='\t',file=Outfile)
main()
