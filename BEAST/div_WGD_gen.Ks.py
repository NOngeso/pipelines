#!/usr/bin/python3 
import numpy,sys
class synteny:
	def __init__(self,strSB):
		self.strSB 		= strSB
		self.Qgenes_list 	= []
		self.Tgenes_list 	= []
		self.Ks_list		= []
		self.fmedianKs		= 0
		self.WGD_pair		= [0,'']
	def parse(self,list_in): #  0-  0:        Glyma01g22260.1 Medtr1g093600.1  2e-142 0.190434613618596       1.61734782410761
		self.Qgenes_list.append(list_in[1])
		self.Tgenes_list.append(list_in[2])
		self.Ks_list.append(float(list_in[5]))
	def calc_medianKs(self):
		self.fmedianKs = numpy.median(self.Ks_list)

def get_div_time(Ks):
	fKs             = float(Ks)
	constant        = 5.17 * 10**(-3)
	divergence_time = fKs/(2*constant)
	return(divergence_time)






file_in = sys.argv[1] #'mt2gm.collinearity.kaks'
Outfile = open(file_in+'.recentpeakWGD','w')
Outfile_target = open(file_in+'.target','w')
targetMin = 0.04 
targetMax = 0.12

dicSB2cls = {}
for line in open(file_in):
	if line[0] == '#':
		continue
	cell = line.strip().split('\t')
	strSB	= cell[0].split('-')[0]
	try:
		cSynt = dicSB2cls[strSB]
	except KeyError:
		cSynt 	= synteny(strSB)
		dicSB2cls[strSB] = cSynt
	cSynt.parse(cell)
dicW2count = {} # window to count
window_start = 0
window_width = 0.01
#border		= float(sys.argv[2]) # 0.05

while int(window_start) != 3:
	window = (window_start,window_start+window_width)
#	print(window)
	dicW2count[window] = 0 
	window_start += window_width

print('#passed')		

target_list = []
for strSB in dicSB2cls:
	cSynt = dicSB2cls[strSB]
	cSynt.calc_medianKs()
	for window in dicW2count:
		if window[0] <= cSynt.fmedianKs < window[1]:
			dicW2count[window] += 1
		if targetMin <= cSynt.fmedianKs <= targetMax:
			target_list.append(cSynt.fmedianKs)
target_list.sort()
print('Ks:',target_list[0],target_list[-1],numpy.average(target_list),file=Outfile_target)
print('MYA:',get_div_time(target_list[0]),get_div_time(target_list[-1]),get_div_time(numpy.average(target_list)),file=Outfile_target)
dicW2count_list = list(dicW2count.keys())
dicW2count_list.sort(key = lambda x : dicW2count[x],reverse=True)
peak_Ks_window = dicW2count_list[0]
print('#','peak Ks window is',peak_Ks_window)
for strSB in dicSB2cls:
	cSynt = dicSB2cls[strSB]
	if strSB != cSynt.strSB:
		print(strSB,cSynt.strSB)
		exit()
	#if peak_Ks_window[0]-border <= cSynt.fmedianKs < peak_Ks_window[1]+border:
	if targetMin <= cSynt.fmedianKs <= targetMax:
		i = 0
		while 1:
			try:
				print(cSynt.strSB,cSynt.Qgenes_list[i],cSynt.Tgenes_list[i],cSynt.Ks_list[i],sep='\t',file=Outfile)
			except IndexError:
				break
			i += 1
dicW2count_list.sort(key = lambda x : x[0])
for window in dicW2count_list:
	print("%0.2f,%0.2f,%d"%(window[0],window[1],int((dicW2count[window]))),'|'* int((dicW2count[window])/10))
	
#	print(cSynt.fmedianKs)


