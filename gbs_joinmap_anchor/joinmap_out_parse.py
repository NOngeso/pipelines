#!/usr/bin/python3

import sys,math

def regression(x,y):

	xsum = 0
	ysum = 0
	xsquaresum = 0
	ysquaresum = 0

	for i in x:
		xsum = xsum + i
		xsquaresum = xsquaresum + i*i

	for j in y:
		ysum = ysum + j
		ysquaresum = ysquaresum + j*j

	xmean = xsum / len(x) # mean of x
	ymean = ysum / len(y) # mean of y

	tsum = 0
	for i in x:
		tsum = tsum + (i-xmean)*(i-xmean) # tsum : temporal variable to sum up numbers
	xd = tsum/(len(x)-1) # variance of x regarding as sampling (N-1)

	tsum = 0
	for j in y:
		tsum = tsum + (j-ymean)*(j-ymean)
	yd = tsum/(len(y)-1) # variance of y regarding as sampling (N-1)

	tsum = 0
	for i in range(0,len(x)):
		tsum = tsum + ((x[i]-xmean) * (y[i]-ymean))

	xd = math.sqrt(xd) # deviation of x
	yd = math.sqrt(yd) # deviation of y
	cv = tsum /(len(y)-1)  # covariance

	residual = cv/(xd*yd) # residual

	return(residual)
	

file_jo = sys.argv[1]#'joinmap_out.txt'
Outfile = open(file_jo+'_ordered.txt','w')
#Outfile_scaf_list = open('joinmap_ordered_scaflist.txt','w')
dicGR2dicSC = {}
for line in open(file_jo):
	if line[0] == ';':
		continue
	if line.strip()== '':
		continue
	if line[0] == '*':
		pass
	elif line[0:5] == 'group':
		cell 	= line.strip().split()
		strGR 	= cell[1]
		dicGR2dicSC[strGR] = {}
		continue
	 
	cell 	= line.strip().split()
	strSC	= ''.join(cell[0].split('_')[0:2])
	intPPOS	= int(cell[0].split('_')[2])
	fCPOS	= float(cell[1])
	try:
		if fCPOS not in dicGR2dicSC[strGR][strSC][1]:
			dicGR2dicSC[strGR][strSC][0].append(intPPOS)
			dicGR2dicSC[strGR][strSC][1].append(fCPOS)
	except KeyError:
		dicGR2dicSC[strGR][strSC] = [[intPPOS],[fCPOS]]

def avr(x):
	tsum = 0
	for i in x:
		tsum += i
	return(tsum/len(x))

for strGR in dicGR2dicSC:
	for strSC in dicGR2dicSC[strGR]:
		x = dicGR2dicSC[strGR][strSC][0]
		y = dicGR2dicSC[strGR][strSC][1]
		
		if len(x) < 2:
			print(strGR,strSC,'N',round(avr(y),2),','.join(map(str,x)),','.join(map(str,y)),sep='\t',file=Outfile)
		else:
			if regression(x,y) > 0:
				strDeter = 'F'
			elif regression(x,y) < 0:
				strDeter = 'R'
			else : strDeter = 'N1'
			print(strGR,strSC,strDeter,round(avr(y),2),','.join(map(str,x)),','.join(map(str,y)),sep='\t',file=Outfile)
	
