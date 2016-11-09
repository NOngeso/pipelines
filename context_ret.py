#!/usr/bin/python
# SNP context retriever 

import sys
sys.path.append('/storage/pipelines/')
import kang


flanklen = 30 

file_fa = sys.argv[1]
dicfa   = kang.Fasta2dic(file_fa)
schr    = sys.argv[2] 
loc     = int(sys.argv[3])-1 ## input location is 1-based position

print dicfa[schr][loc-flanklen:loc],dicfa[schr][loc],dicfa[schr][loc+1:loc+1+flanklen]

