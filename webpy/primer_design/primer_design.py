#!/usr/bin/python

import web  
import sys
sys.path.append('../') 
import kang

file_fa = 'Creinhardtii_281_v5.0.fa'
dicHD2seq = kang.Fasta2dic(file_fa)

urls = ('/', 'index') # . ... .. (/) . .... index.. .... ... ....

render = web.template.render('templates/') # ...... .... html... ... ... .... 

