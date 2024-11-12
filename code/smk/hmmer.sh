# Description: run hmmsearch vs all proteins
for i in /data/san/data1/users/david/mining_2023/models/nis*.hmm; do
    outfile=$(basename $i) 
	hmmsearch --cpu 8 --tblout /data/san/data2/users/david/nisin/data/tables/hmm/${outfile%.hmm}_hmm.out $i /data/san/data2/users/david/nisin/data/all_proteins_ncbi/all_proteins2.faa
done
