Q=Athaliana_167_protein.fa
T=Alyrata_107_protein.fa
N=At2Al
blastall -i ${Q} -d ${T} -p blastp -e 1e-10 -b 5 -v 5 -m8 -o ${N}.blast.pre -a 5
