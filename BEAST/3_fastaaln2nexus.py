#!/usr/bin/python3
'''
#NEXUS
Begin data;
Dimensions ntax=4 nchar=15;
Format datatype=dna symbols="ACTG" missing=? gap=-;
Matrix
Species1   atgctagctagctcg
Species2   atgcta??tag-tag
Species3   atgttagctag-tgg
Species4   atgttagctag-tag           
;
End;
'''
import kang,glob

dicAcc2tax = {'VMM':'Ceratotropis','VMS':'Ceratotropis','VGR':'Ceratotropis','VRA':'Ceratotropis','VRS':'Ceratotropis','VSU':'Ceratotropis','VAC':'Aconitifoliae','VTR':'Aconitifoliae','VST':'Aconitifoliae','VHA':'Aconitifoliae','VUC':'Angulares','VU2':'Angulares','VNE':'Angulares','VAG':'Angulares','VNA':'Angulares','VRI':'Angulares','VMI':'Angulares','VRG':'Angulares','VRR':'Angulares','VTA':'Angulares','V38':'Eurasian','V37':'African'}

def fasta2nexus(dicHD2seq,formattype): 
	dicHD2seq_list = list(dicHD2seq.keys())
	dicF2sym = {}
	dicF2sym['protein'] 	= '"FSTNKEYVQMCLAWPHDRIG"'
	dicF2sym['dna']		= '"ACTG"'
	result = '''#NEXUS
Begin data;
Dimensions ntax=%d nchar=%d;
Format datatype=%s symbols=%s missing=? gap=-;
Matrix
'''%(len(dicHD2seq_list),len(dicHD2seq[dicHD2seq_list[0]]),formattype,dicF2sym[formattype])
	for key in dicHD2seq:
	#	try:
#			acc 	= dicAcc2tax[key.split('|')[0]] + '|' + key  # VRG|m.86413_1
#		except KeyError:
#			acc     = key + '|' + key.split('|')[0] # VRG|m.86413_1
		acc	= key
		seq 	= dicHD2seq[key]
		result += acc + '\t' + seq + '\n'
	result += ';\n'
	result += 'End;\n'
	return(result)
fas_list	= glob.glob('*.best.fas')
for file_fas in fas_list:
	dicHD2seq 	= kang.Fasta2dic(file_fas)
	nexus 		= fasta2nexus(dicHD2seq,'protein')
	Outfile = open('.'.join(file_fas.split('.')[0:])+'.nex','w')
	print(nexus.strip(),file=Outfile)
	





	
