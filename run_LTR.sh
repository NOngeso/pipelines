/data2/k821209/programs/genometools-1.5.1/bin/gt suffixerator -suf -lcp -dna -indexname ${1} -db ${1}
/data2/k821209/programs/genometools-1.5.1/bin/gt ltrharvest -gff3 ${1}.gff3 -index ${1}
/data2/k821209/programs/genometools-1.5.1/bin/gt gff3 -sort ${1}.gff3 > ${1}.sort.gff3
/data2/k821209/programs/genometools-1.5.1/bin/gt -j 2 ltrdigest -pptlen 10 30 -pbsoffset 0 3 -trnas Athal-tRNAs.fa -hmms ./hmm/*.hmm -outfileprefix ${1} ${1}.sort.gff3 Athaliana_167.fa | tee ${1}.ltrdigest_GyDB.gff3

