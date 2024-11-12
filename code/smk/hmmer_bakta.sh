# write a loop to run over all hmm files and do hmmsearch vs all proteins

# Path: hmmer.sh
# Author: <NAME>
# Date: 2019-01-29
# Description: run hmmsearch vs all proteins

for i in /data/san/data1/users/david/mining_2023/models/nis*.hmm; do
    outfile=$(basename $i) # This gets just the filename, without the path
	hmmsearch --cpu 8 --tblout /data/san/data2/users/david/nisin/data/tables/hmm_bakta/${outfile%.hmm}_hmm.out $i /data/san/data2/users/david/nisin/data/ava/all_proteins.faa
done
