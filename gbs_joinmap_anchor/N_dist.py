#!/usr/bin/python3


intNN_cutoff	= 30
merge_range	= 1000
file_sc 	= 'total.txt'
Outfile_dist 	= open(file_sc+'.dist','w')
Outfile		= open(file_sc+'.%dfilt.merge%d.represent'%(intNN_cutoff,merge_range),'w')
dicNN2NS = {} # Number of 'N' to Number of SNP locations
def genotyping(strREF,GT_list):
	stack = ''
	for each in GT_list:
		if each == strREF:
			stack += 'a'
		elif each == 'N':
			stack += '-'
		elif each == 'H':
			stack += 'h'
		else : stack += 'b'
	return(stack)
haplotype = []

def merge(haplotype):
	genotype_list = [x[2] for x in haplotype]
	merged = ''
	for i in range(len(genotype_list[0])):
		stack = []
		for genotype in genotype_list:
			stack.append(genotype[i])
		'''
		if 'a' in stack and 'b' in stack:
			merged += '-'
		elif 'a' in stack and 'b' not in stack:
			merged += 'a'
		elif 'a' not in stack and 'b' in stack:
			merged += 'b'
		elif 'a' not in stack and 'b' not in stack:
			merged += '-'
		'''
		iA = stack.count('a')
		iB = stack.count('b')
		iH = stack.count('h')
		iN = stack.count('-')
		tot = iA + iB + iH + iN
		'''
		if iA == 0 and iB == 0 and iH > 0:
			merged += 'h'
		elif iA == 0 and iB > 0 and iH == 0:
			merged += 'b'
		elif iA > 0 and iB == 0 and iH == 0:
			merged += 'a'
		else: merged += '-'
		'''
#		print(stack)
		if iA / tot > 0.8:
			merged += 'a'
		elif iB / tot > 0.8:
			merged += 'b'
		elif iH / tot > 0.8:
			merged += 'h'
		else: merged += '-'
	return(merged)
		
		
		

for line in open(file_sc):
	cell 		= line.strip().split('\t')
	strSC		= cell[0].split('_')[0] + '_' + cell[0].split('_')[1] 
	strLoc		= str(int(cell[0].split('_')[2]) + int(cell[1]))
	strREF		= cell[2]
	strVAR_list	= cell[3].split(',')
	intNN		= strVAR_list.count('N')
	if intNN <= intNN_cutoff:
		#print(strSC,strLoc,strREF,','.join(strVAR_list),file=Outfile,sep='\t')
		if haplotype == []:
			haplotype = [[strSC,strLoc,genotyping(strREF,strVAR_list),intNN]]
		elif haplotype[-1][0] == strSC and abs(int(haplotype[-1][1])-int(strLoc)) < merge_range: # Haplotype is 10k
			haplotype.append([strSC,strLoc,genotyping(strREF,strVAR_list),intNN])
		else : 
			haplotype.sort(key = lambda x:x[3])
			genotypes = merge(haplotype)
			if genotypes.count('-') < intNN_cutoff: # when surrounded genotype is not consistent each other, skip
				print('*'+haplotype[0][0].replace('caffold','').replace('uperScaf','S') + '_' + haplotype[0][1]+'_'+str(len(haplotype)),'\t'.join(list(genotypes)),file=Outfile,sep='\t')
			haplotype = [[strSC,strLoc,genotyping(strREF,strVAR_list),intNN]]
		#print(strSC+'_'+strLoc,genotyping(strREF,strVAR_list),file=Outfile)
	try:
		dicNN2NS[intNN] += 1
	except KeyError:
		dicNN2NS[intNN] = 1
dicNN2NS_list = list(dicNN2NS.keys())
dicNN2NS_list.sort()
stack = 0
for intNN in dicNN2NS_list:
	stack += dicNN2NS[intNN]
	print(intNN,stack,file=Outfile_dist)

	
	
