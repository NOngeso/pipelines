#!/usr/bin/python3

import sys
#########################################################################
#strEnzyme_site_list = [['PstI','CTGCAG',5],['ApeKI','GCWGC',1],['XbaI','TCTAGA',1]]
#strEnzyme_site_list = [['HindIII','AAGCTT',1],['NlaIII','CATG',4]]
strEnzyme_site_list = [['ApeKI','GCWGC',1]]
Outfile_block       = open('Enzyme_site_block_ApeKI.txt','w')
#########################################################################
file_fasta = sys.argv[1] #'../final.assembly.fasta'
#file_fasta = 'test.fa'
bulk_load  = open(file_fasta).read()
bulk_list  = bulk_load.split('>')
dicHD2Seq  = {}
for each_bulk in bulk_list:
        if each_bulk == '':
                continue
        strHD            = each_bulk.split('\n')[0]
        strSeq           = ''.join(each_bulk.split('\n')[1:])
        dicHD2Seq[strHD] = strSeq
#########################################################################


def check_enz(strChr,intPos):
	for strEz_name, strEz_seq, intCutPos in strEnzyme_site_list:
		if intPos+len(strEz_seq) > len(dicHD2Seq[strChr]):
			continue
		if 'W' in strEz_seq:
			strEz_seq1 = strEz_seq.replace('W','A')
			strEz_seq2 = strEz_seq.replace('W','T')
			if dicHD2Seq[strChr][intPos:intPos+len(strEz_seq)] == strEz_seq1 or dicHD2Seq[strChr][intPos:intPos+len(strEz_seq)] == strEz_seq2:
				return(1,intCutPos,strEz_name,strEz_seq)
		elif dicHD2Seq[strChr][intPos:intPos+len(strEz_seq)] == strEz_seq:
			return(1,intCutPos,strEz_name,strEz_seq)
	return(0,0,0,0)
def main():
	intTotal_subgenomesize = 0
	for dicHD2Seq_key in dicHD2Seq:
		dicHD2Seq_value = dicHD2Seq[dicHD2Seq_key]
		Pre_pos  = -1
		Pre_name = ''
		Pre_seq  = ''
		for i in range(len(dicHD2Seq_value)):
			bExist, intPos, strEz_name, strEz_seq = check_enz(dicHD2Seq_key,i)
			if bExist == 1:
				intCutPos = i + intPos
				if Pre_pos != '-1':
					if 2000 > len(dicHD2Seq_value[Pre_pos:intCutPos]) > 10:
						if Pre_name == 'ApeKI' and strEz_name == 'ApeKI':
							#print(dicHD2Seq_key, Pre_pos, intCutPos, Pre_name+':'+'['+Pre_seq+']','~',strEz_name+':'+'['+strEz_seq+']',len(dicHD2Seq_value[Pre_pos:intCutPos]))
							print(dicHD2Seq_key, Pre_pos, intCutPos, Pre_name+':'+'['+Pre_seq+']','~',strEz_name+':'+'['+strEz_seq+']',len(dicHD2Seq_value[Pre_pos:intCutPos]),dicHD2Seq_value[Pre_pos-5:Pre_pos],dicHD2Seq_value[Pre_pos:intCutPos],dicHD2Seq_value[intCutPos:intCutPos+5],file=Outfile_block)
							intTotal_subgenomesize += len(dicHD2Seq_value[Pre_pos:intCutPos])
					Pre_pos  = intCutPos
					Pre_name = strEz_name
					Pre_seq  = strEz_seq
				if Pre_pos == '-1':
					Pre_pos  = intCutPos
					Pre_name = strEz_name
					Pre_seq  = strEz_seq
	print('## Total subgenome size = %d'%intTotal_subgenomesize,file = Outfile_block)
main()


	
