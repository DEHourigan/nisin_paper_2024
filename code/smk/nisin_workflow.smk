# write a snakfile to run gtdbtk on a folder of genomes tha
# is in /data/san/data1/users/david/mining_2023/nisins/genomes_download/assemblies
# and output the results to /data/san/data1/users/david/mining_2023/nisins/genomes_download/gtdb

# the script should be able to run on the cluster with the following command
# the input is a single folder containing all genomes
# out put is a series of directories containing the gtdbtk 
import pandas as pd
import argparse
import glob

def get_samples():
	return [x.split('/')[-1].replace('.fna', '') for x in glob.glob("../../data/bakta_out/*/*.fna")]


rule all:
	input:
		# "gtdb/classify/classify/gtdbtk.bac120.summary.tsv",
		'/data/san/data2/users/david/nisin/data/tables/genbank_meta.tsv',
		'/data/san/data2/users/david/nisin/data/tables/assemblies_genbank.txt',
		'/data/san/data2/users/david/nisin/data/tables/geo_loc_name.tsv',
		'/data/san/data2/users/david/nisin/data/tables/isolation_source.tsv',
		'/data/san/data2/users/david/nisin/data/tables/host.tsv',
		'/data/san/data2/users/david/nisin/data/tables/collection_date.tsv',
		'/data/san/data2/users/david/nisin/data/tables/hmm/nisT_hmm.out',
		expand("../../data/repeat/{sample}.2.5.7.80.10.50.2000.summary.html", sample=get_samples()),
		# expand("../../data/repeat/{sample}.s2.5.7.80.10.50.2000.txt.html", sample=get_samples()),
		# expand("../../data/repeat/{sample}.s2.5.7.80.10.50.2000.html", sample=get_samples()),
		# expand("../../data/repeat/{sample}.s1.5.7.80.10.50.2000.txt.html", sample=get_samples()),
		# expand("../../data/repeat/{sample}.s1.5.7.80.10.50.2000.html", sample=get_samples()),
		expand('../../data/tncomp/{sample}/blastn/{sample}.blastn', sample=get_samples())



		
# rule gtdbtk:
# 	workdir: "/data/san/data1/users/david/mining_2023/nisins/genomes_download"
# 	input:
# 		"assemblies/"v
# 	output:
# 		"gtdb"
# 	conda:
# 		"/home/david/yamls/gtdbtk23.yml"
# 	shell:
# 		"gtdbtk identify --genome_dir {input} --out_dir {output}/identify"
# 		"gtdbtk align --identify_dir {output}/identify --out_dir {output}/align"
# 		"gtdbtk classify --identify_dir {output}/identify --out_dir {output}/classify"
# 		"gtdbtk infer --identify_dir {output}/identify --out_dir {output}/infer"

rule mv_acc:
	input:
		"/data/san/data1/users/david/mining_2023/nisins/genomes_download/assemblies_genbank.txt"
	output:
		"/data/san/data2/users/david/nisin/data/tables/assemblies_genbank.txt"
	shell:
		"""
		cp {input} {output}
		"""

# use datasets to get meta data for the Genbank assemblies
rule get_meta:
	input:
		"/data/san/data2/users/david/nisin/data/tables/assemblies_genbank.txt"
	output:
		"/data/san/data2/users/david/nisin/data/tables/genbank_meta.tsv"
	conda: "/home/david/yamls/datasets.yml"
	shell:
		"""
		datasets summary genome accession --inputfile {input} --as-json-lines | dataformat tsv genome --fields accession,assminfo-biosample-accession,assmstats-gc-percent,checkm-version,checkm-completeness,assminfo-biosample-description-organism-infraspecific-isolate,assmstats-total-sequence-len,assminfo-biosample-attribute-name,assminfo-biosample-attribute-value > {output}
		"""

# filter the table genbank_meta.tsv to make a table called isolation country
rule filter_attributes:
	input:
		infile="/data/san/data2/users/david/nisin/data/tables/genbank_meta.tsv"
	output:
		country='/data/san/data2/users/david/nisin/data/tables/geo_loc_name.tsv',
		source="/data/san/data2/users/david/nisin/data/tables/isolation_source.tsv",
		host="/data/san/data2/users/david/nisin/data/tables/host.tsv",
		date="/data/san/data2/users/david/nisin/data/tables/collection_date.tsv"
	params:
		country="geo_loc_name",  # Change this to the desired attribute name
		source="isolation_source",
		host="host",
		date="collection_date"

	shell:
		"""
			python ../meta_table_sort.py --infile {input.infile} --param {params.country} --outfile {output.country}
			python ../meta_table_sort.py --infile {input.infile} --param {params.source} --outfile {output.source}
			python ../meta_table_sort.py --infile {input.infile} --param {params.host} --outfile {output.host}
			python ../meta_table_sort.py --infile {input.infile} --param {params.date} --outfile {output.date}			
		"""

# FILEPATH: /data/san/data2/users/david/nisin/code/smk/nisin_workflow.smk

rule repeats:
	input:
		fna_file = "../../data/bakta_out/{sample}/{sample}.fna"
	output:
		out = "../../data/repeat/{sample}.2.5.7.80.10.50.2000.summary.html",
		outdir="../../data/repeat/"
	threads:
		12
	shell:
		"""
		source ~/.profile 
		mkdir -p ../../data/repeat
		echo "Running Repeat Finder on {input.fna_file}"
		cd {output.outdir}
		trf {input.fna_file} 2 5 7 80 10 50 2000 -l 1 -ngs > {output.out}
		"""

# rule icefinder:
#     input:
#         fna_file = "../../data/bakta_out/{sample}/{sample}.fna"
#     output:
#         dummy = "../../data/ice/ICEfinder2_linux/result/{sample}.dummy"
#     threads:
#         12
#     conda:
#         "/home/david/yamls/icefinder.yml"
#     shell:
#         """
#         echo "Running IceFinder on {input.fna_file}"
#         python ../../data/ice/ICEfinder2_linux/ICEfinder2.py -i {input.fna_file} -t Single
#         touch {output.dummy}
#         """

rule tncomp:
    input:
        fna_file = "../../data/tncomp/{sample}.fna"
    params: 
        '../../data/tncomp/{sample}_out'
    threads:
        12
    conda:
        "/home/david/yamls/tncomp.yml"
    shell:
        """
        python /data/san/data2/users/david/programmes/tncomp_finder/TnComp_finder.py -f {input.fna_file} -o {params} -p 16
        """


# rule plasmid_content:
# 	input:
# 		fna_file = "../../data/bakta_out/{sample}/{sample}.fna"
# 	output:
# 		out_dir = "../../data/tncomp/{sample}/"
# 	threads:
# 		12
# 	conda:
# 		"/home/david/yamls/tncomp.yml"
# 	shell:
# 		"""
# 		echo "Running tncomp on {wildcards.sample}"
# 		mkdir -p {output.out_dir}
# 		python /data/san/data2/users/david/programmes/tncomp_finder/TnComp_finder.py -f {input.fna_file} -o {output.out_dir} -p 16
# 		"""

