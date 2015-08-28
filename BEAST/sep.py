#!/usr/bin/python3 

import numpy,kang

file_fa	= 'final.assembly.11.fasta.nonATGC.fa.cds.primary.fasta.pep.sh.fa'
file_vra2vrp = 'Vra2Vrr.collinearity.kaks.recentpeakWGD'
file_vrp2vrp = 'Vrg2Vrg.collinearity.kaks.recentpeakWGD'

dicHD2seq = kang.Fasta2dic(file_fa)

dicGN2info = {}

for line in open(file_vra2vrp):
	if line[0] == '#':
		continue
	cell = line.strip().split('\t')
	strSB	= cell[0].split('-')[0]
	strGN	= cell[2]
	strKs	= cell[3]
	try:
		if dicGN2info[strGN]:
			print(strGN)
	except KeyError:
		dicGN2info[strGN] = [strSB,strKs]

class sb():
	def __init__(self):
		self.strSB = ''
		self.G1_Ks_list = []
		self.G1_SB_list = []
		self.G1_list	= []
		self.G2_Ks_list = []
		self.G2_SB_list = []
		self.G2_list	= []
		
dicSB2c = {}
for line in open(file_vrp2vrp):
	cell 	= line.strip().split('\t')
	strSB 	= cell[0]
	strG1	= cell[1]
	strG2	= cell[2]
	strKs	= cell[3]
	try:
		cSb = dicSB2c[strSB]
		try:
			cSB.G1_Ks_list.append(float(dicGN2info[strG1][1]))
			cSB.G1_SB_list.append(dicGN2info[strG1][0])
			cSB.G2_Ks_list.append(float(dicGN2info[strG2][1]))
			cSB.G2_SB_list.append(dicGN2info[strG2][0])
			cSB.G1_list.append(strG1)
			cSB.G2_list.append(strG2)
		except KeyError:
			pass
	except KeyError:
		cSB 	= sb()
		dicSB2c[strSB] = cSB
		cSB.strSB = strSB
		try:
			cSB.G1_Ks_list.append(float(dicGN2info[strG1][1]))
			cSB.G1_SB_list.append(dicGN2info[strG1][0])
			cSB.G2_Ks_list.append(float(dicGN2info[strG2][1]))
			cSB.G2_SB_list.append(dicGN2info[strG2][0])
			cSB.G1_list.append(strG1)
			cSB.G2_list.append(strG2)
		except KeyError:
			pass
print(1)
VRPA = []
VRPB = []
VRPA_ks = []
VRPB_ks = []	
for strSB in dicSB2c:
	cSB = dicSB2c[strSB]
	if len(set(cSB.G1_SB_list)) != 1 or len(set(cSB.G2_SB_list)) != 1:
		#rint(cSB.G1_SB_list)
		continue
	#print(numpy.median(cSB.G1_Ks_list),numpy.median(cSB.G2_Ks_list))
	G1_mKs = numpy.median(cSB.G1_Ks_list)
	G2_mKs = numpy.median(cSB.G2_Ks_list)
	if G1_mKs < G2_mKs:
		VRPA += cSB.G1_list
		VRPB += cSB.G2_list
		VRPA_ks += cSB.G1_Ks_list
		VRPB_ks += cSB.G2_Ks_list
	else:
		VRPA += cSB.G2_list
		VRPB += cSB.G1_list
		VRPA_ks += cSB.G2_Ks_list
		VRPB_ks += cSB.G1_Ks_list
#print(VRPA)
#print(VRPB)
print(numpy.median(VRPA_ks),numpy.median(VRPB_ks))
Outfile = open('VRPA.fasta','w')
for gene in VRPA:
	print('>'+gene,file=Outfile)
	print(dicHD2seq[gene],file=Outfile)
Outfile.close()
Outfile = open('VRPB.fasta','w')
for gene in VRPB:
	print('>'+gene,file=Outfile)
	print(dicHD2seq[gene],file=Outfile)
Outfile.close()
Outfile = open('match.txt','w')
for i in range(len(VRPA)):
	print(VRPA[i],VRPB[i],sep='\t',file=Outfile)
Outfile.close()
