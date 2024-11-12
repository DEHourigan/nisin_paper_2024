import glob
import os

# Get all .fna files in the input directory
fna_files = glob.glob("../../data/tncomp/in/*.fna")
# Extract sample names without extensions and paths
samples = [os.path.basename(f).replace(".fna", "") for f in fna_files]

# Print debug information
print("Found samples:", samples)

# Define the rule all to specify the final output targets
rule all:
    input:
        expand("../../data/tncomp/{sample}_out", sample=samples)

# Define the tncomp rule
rule tncomp:
    input:
        fna_file = "../../data/tncomp/in/{sample}.fna"
    params:
        out_dir = "../../data/tncomp/{sample}_out"
    threads:
        12
    conda:
        "/home/david/yamls/tncomp.yml"
    shell:
        """
        echo "Running TnComp on {input.fna_file} with output directory {params.out_dir}"
        python /data/san/data2/users/david/programmes/tncomp_finder/TnComp_finder.py -f {input.fna_file} -o {params.out_dir} -p 16
        """
