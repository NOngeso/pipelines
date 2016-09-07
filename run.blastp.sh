Q=${1}
S=${2}
blastp -query $Q -db $S -evalue 1e-5 -num_alignments 1 -outfmt 7 -out $Q.bp.ev1e5.na1.$S.out7 -num_threads 5
