from __future__ import print_function
import subprocess
import numpy as np
import pandas as pd
import sys
sys.path.append('/ref/analysis/pipelines/')

import kang
file_bam = 'merge.all.bam'
file_fa = '../../ref/Creinhardtii_281_v5.0.fa'
dicHD2seq = kang.Fasta2dic(file_fa)

#chromosome, left, right = 'chromosome_1', 1000, 20000

#print ('samtools view %s %s:%d-%d > temp.sam'%(file_bam,chromosome,left,right))
#subprocess.call('samtools view %s %s:%d-%d > temp.sam'%(file_bam,chromosome,left,right),shell=True)
#subprocess.call('cut -f 1-20 temp.sam > temp.sam.cut',shell=True)


rows              = len(dicHD2seq.keys())
chromosomes       = dicHD2seq.keys()
chromosomes.sort()
dicN2chr          = dict(enumerate(chromosomes))
dicChr2N          = {b:a for a,b in dicN2chr.iteritems()}
columns           = max([len(x) for x in dicHD2seq.values()])-1
continuity_matrix = np.zeros([rows,columns],dtype=np.int)
Outfile = open('chromosome.map.txt','w')
for a,b in dicChr2N.iteritems():
    print(a,b,sep='\t',file=Outfile)


print('start loop')
i = 0 
for line in sys.stdin:#$open('temp.sam.cut'): # should be changed to zero base map
    #print line
    cell         = line.strip().split('\t')
    chromosome   = cell[2]
    echr         = dicChr2N[chromosome]
    startpos     = int(cell[3]) - 1 # zero based conversion 
    fragmentsize = int(cell[8])
    if fragmentsize > 200:
        pass
    else: continue
    if i % 100000 == 0 :
        print(i, 'th line processing')
    i += 1
    endpos       = startpos + fragmentsize - 1 
    startpos_inter = startpos  # startpos - nextpos 
    endpos_inter   = endpos - 1 # previouspos - endpos 
    
    continuity_matrix[echr,startpos_inter:fragmentsize-1] += 1  # list characteristic can utillize fragment size itself.
    
print('saving..')
np.save('test.np',continuity_matrix)

#np.savetxt('test.txt',continuity_matrix)    


