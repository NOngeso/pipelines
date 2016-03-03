from __future__ import print_function


gencode = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'}
def rev_comp(strSeq):
	dicComp = {'A':'T','T':'A','G':'C','C':'G','N':'N'}
	strCseq = ''
	for i in strSeq:
		try:
			strCseq += dicComp[i.upper()]
		except KeyError:
			strCseq += 'N'
	# End of for i
	return(strCseq[::-1])
def translation(strSeq):
	strPep = ''
	for i in range(int(len(strSeq)/3)):
		try:
			strPep += gencode[strSeq[i*3:i*3+3].upper()]
		except KeyError:
			strPep += 'X'
	# End of for i
	return(strPep)
def fasta2dic(file_fasta, dic):
	bulk      = open(file_fasta).read()
	bulk_list = bulk.split('>')
	for each in bulk_list:
		if each == '':
			continue
		strHD  = each.split('\n')[0].split()[0]
		strSeq = ''.join(each.split('\n')[1:])
		dic[strHD] = strSeq
	return(dic)
def Fasta2dic(file_fasta):
	dic       = {}
	bulk      = open(file_fasta).read()
	bulk_list = bulk.split('>')
	for each in bulk_list:
		if each == '':
			continue
		strHD  = each.split('\n')[0].split()[0]
		strSeq = ''.join(each.split('\n')[1:])
		dic[strHD] = strSeq
	return(dic)
def Fasta2dic_all(file_fasta):
	dic       = {}
	bulk      = open(file_fasta).read()
	bulk_list = bulk.split('>')
	for each in bulk_list:
		if each == '':
			continue
		strHD  = each.split('\n')[0]
		strSeq = ''.join(each.split('\n')[1:])
		dic[strHD] = strSeq
	return(dic)

def dic2fa(dic,filename):
	Outfile = open(filename,'w')
	for key in dic:
		print('>'+key,file=Outfile)
		print(dic[key],file=Outfile)
	Outfile.close()

		



