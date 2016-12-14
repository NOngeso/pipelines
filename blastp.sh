Q=Athaliana_167_TAIR10.protein.fa
D=Athaliana_167_TAIR10.protein.fa
blastp -query ${Q} -db ${D} -evalue 1e-5 -num_threads 5 -outfmt 7 -out ${Q}.${D}.bp.ev1e5.out7 
