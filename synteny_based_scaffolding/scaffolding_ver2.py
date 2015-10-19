#!/usr/bin/python3

import kang

class SB():
	def __init__(self,strSB):
		self.strSB = strSB
		self.strSC_order = []
		self.dicSC2LOC = {}
Outfile = open('scaffolds.txt','w')
Outfile_trace = open('scaffolds.trace','w')


SubSB 	= 'Pv2Vr.collinearity.kaks.recentpeakWGD'
MainSB	= 'Pv2Mt.collinearity.kaks.recentpeakWGD'
Vagff	= 'vr.gff'
dicVA2SC	= {}
dicVA2LOC 	= {}
for line in open(Vagff):
	cell 	= line.strip().split('\t')
	strSC 	= cell[0]
	strGN	= cell[1]
	strLOC	= cell[2]
	dicVA2SC[strGN] 	= strSC
	dicVA2LOC[strGN] 	= strLOC

done = []
dicMainSB2MGN	= {}
for line in open(MainSB):
	cell 	= line.strip().split('\t')
	strSB	= cell[0] 
	strMGN	= cell[1] # MGN : Medium genename 
	try:
		dicMainSB2MGN[strSB].append(strMGN)
	except KeyError:
		dicMainSB2MGN[strSB] = [strMGN]

dicMGN2SubSB = {}
dicSubSB2info = {}
for line in open(SubSB):
	cell	= line.strip().split('\t')
	strSB	= cell[0]
	strMGN	= cell[1]
	strGN	= cell[2]
	try:
		dicSubSB2info[strSB][0].append(dicVA2SC[strGN])
		dicSubSB2info[strSB][1].append(dicVA2LOC[strGN])
	except KeyError:
		dicSubSB2info[strSB] = [[dicVA2SC[strGN]],[dicVA2LOC[strGN]]]
	try:
		dicMGN2SubSB[strMGN][0].append(strSB)
		dicMGN2SubSB[strMGN][1].append(strGN)
	except KeyError:
		dicMGN2SubSB[strMGN] = [[strSB],[strGN]]

def stackzip(list1,list2): # length of two lists should be identical
	dic = {}
	for i in range(len(list1)):
		try:
			dic[list1[i]].append(list2[i])
		except KeyError:
			dic[list1[i]] = [list2[i]]
	return(dic)

dicSB2cSB = {}

for strSB in dicMainSB2MGN:
	subSB_list 	= []
	dicSubSB2num 	= {}
	for strMGN in dicMainSB2MGN[strSB]:
		try:
			subSBs = dicMGN2SubSB[strMGN][0]
			subGenes = dicMGN2SubSB[strMGN][1]
		except KeyError:
			continue
		for subSB in subSBs:
			try:
				dicSubSB2num[subSB] += 1
			except KeyError:
				dicSubSB2num[subSB] = 1
	subSB_list = [x for x in dicSubSB2num if dicSubSB2num[x] > 15] # important : number of genes to determine multi-species synteny relationship
	# subSB_list represent the matehced Va2Vr synteny blocks for the Vr2GM syntent block
	# then, how the sub synteny blocks were mapped to Vr genome? 
	# The order of the synteny should be marked.
	if len(subSB_list) > 1:
		#print(strSB,SB_list,dicSB2VrGm[strSB])
		try:
			cSB = dicSB2cSB[strSB]
		except KeyError:
			cSB = SB(strSB)
			dicSB2cSB[strSB] = cSB
		#cSB.strSC_order = [dicSubSB2info[x][0] for x in subSB_list]
		SC_list 	= []
		LOC_list 	= []
		for subSB in subSB_list:
			SC_list		+= dicSubSB2info[subSB][0]
			LOC_list 	+= dicSubSB2info[subSB][1]
		cSB.strSC_order = list(set(SC_list))
		cSB.dicSC2LOC 		= stackzip(SC_list,LOC_list)
def deter_increasing(list_in):
        bINCREASE = []
        for i in range(len(list_in)):
                try:
                        if list_in[i] < list_in[i+1]:
                                bINCREASE.append(1)
                        elif list_in[i] > list_in[i+1]:
                                bINCREASE.append(0)
                except IndexError:
                        pass
        print(list_in)
        print(bINCREASE)
        if sum(bINCREASE) == len(list_in)-1:
		#return('Forward(%d)'%len(list_in)) # Increasing
                return('F')
        elif sum(bINCREASE) == 0:
                #return('Reverse(%d)'%len(list_in)) # Decreasing
                return('R')
        else:
                #print(list_in)
                return('W') # Weird



def rev(strSTRD):
	if strSTRD == 'F':
		return('R')
	elif strSTRD == 'R':
		return('F')
	else:
		return('W')

dicL1toL2 	= {}
edges		= []
vertices 	= []
for strSB in dicSB2cSB:

	cSB = dicSB2cSB[strSB]
	strSC_list = []
	for strSC in cSB.strSC_order:
		print(strSB,strSC,cSB.dicSC2LOC[strSC],deter_increasing(list(map(int,cSB.dicSC2LOC[strSC]))))
		strSC_list.append(strSC+','+deter_increasing(list(map(int,cSB.dicSC2LOC[strSC]))))
	if len(strSC_list) > 1:
		pass
	else : continue
	for i in range(len(strSC_list)):	
		L1 = strSC_list[i]
		try:
			L2 = strSC_list[i+1]
		except IndexError:
			continue
	#	print(L1,L2)
		if [L1,L2] not in edges:
			edges.append([L1,L2])
		if L1 not in vertices:
			vertices.append(L1.split(',')[0])
		if L2 not in vertices:
			vertices.append(L2.split(',')[0])
		try:
			if L2 not in dicL1toL2[L1]:
				dicL1toL2[L1].append(L2)
		except KeyError:
			dicL1toL2[L1] = [L2]
		#Reverse 
		Lr1 = ','.join([L2.split(',')[0],rev(L2.split(',')[1])])
		Lr2 = ','.join([L1.split(',')[0],rev(L1.split(',')[1])])
		#print(Lr1,Lr2)
		if [Lr1,Lr2] not in edges:
			edges.append([Lr1,Lr2])
		try:
			if Lr2 not in dicL1toL2[Lr1]:
				dicL1toL2[Lr1].append(Lr2)
		except KeyError:
			dicL1toL2[Lr1] = [Lr2]
#for L1 in dicL1toL2:
	#print(L1, dicL1toL2[L1])
print(len(vertices))
print(len(edges))
links = []
for edge in edges:
	V1, V2 = edge
	done_edge = []
	stack = []
	stack_list = []
	while stack != [edge[0],edge[1]]: # until no more extension
		pre_len = 0
		V1, V2 = edge
		stack = [V1,V2]
		last_tedge = []
		done_scaf = [V1.split(',')[0],V2.split(',')[0]]
		while 1:
			if len(stack) == pre_len:# last edge remove	
				done_edge.append(last_tedge)
				break
		#	print(len(stack),pre_len)
		#	print(stack)
			pre_len = len(stack)
			for tedge in edges:
				if tedge in done_edge:
		#			print('#',tedge,'#')
					continue
				tV1, tV2 = tedge
				if V2 == tV1:
					if tV2.split(',')[0] not in done_scaf:
						stack += [tV2]
						done_scaf.append(tV2.split(',')[0])
						V2 = tV2
						last_tedge = tedge
						#done_edge.append(tedge)
						break
					else: 
						#rint(tV1,tV2)
						continue
		print(stack,file=Outfile_trace)
		stack_list.append(stack)
		
	print('##',file=Outfile_trace)
	stack_list.sort(key=lambda x:len(x),reverse=True)
#	print(stack_list[0])
	links.append(stack_list[0])

done_ver = []
for vertex in vertices:
	if vertex in done_ver:
		continue
	vertex_links = []
	for link in links:
		if vertex in [x.split(',')[0] for x in link]:
			vertex_links.append(link)
	vertex_links.sort(key=lambda x:len(x),reverse=True)
	member_links = []
	for each in vertex_links[0]:
		for link in links:
			if each in link:
				member_links.append(link)
	member_links.sort(key=lambda x:len(x),reverse=True)	
	bGO = 1
	for each in member_links[0]:
		if each.split(',')[0] in done_ver:
			bGO = 0
	if bGO == 0:
		continue
	print('\t'.join(member_links[0]),file=Outfile)
	for evertex in [x.split(',')[0] for x in member_links[0]]:
		if evertex not in done_ver:
			done_ver.append(evertex)
	
