##
import gzip
import os

def sequencing_bam_file(wildcards):
    individuals = get_individuals()
    sequencing_bam_file_per_individual = individuals[wildcards.individual] # This step here is crucial! It links the individual ID with the bam filename from the individuals.txt file. The latter doesn't even have the individual ID in its file name!

    return expand("{sequencing_bam_dir}/{sequencing_bam_file_per_individual}", sequencing_bam_dir = config["sequencing_bam_dir"], sequencing_bam_file_per_individual = sequencing_bam_file_per_individual )

def sorted_bam_files(wildcards):
    return expand("{bam_dir}/{individual}.ML_included.sorted.bam", bam_dir = config["bam_dir"], individual = wildcards.individual)

rule bams:
    input:
        sorted_bam_files

rule index_reference:
    input: 
        config["genome"]
    output:
        "results/genome/genome.mmi"
    log: expand("{logs}/index_reference.log", logs=config["log_dir"])
    shell:
        "/user/hote1586/tools/minimap2/minimap2 -t 8 -k 17 -I 8G -d {output} {input[0]} > {log} 2>&1"

rule align:
    input:
        sequencing_bam_file,
        "results/genome/genome.mmi"
    output:
        expand("{bam_dir}/{{individual}}.ML_included.sorted.bam", bam_dir = config["bam_dir"])
    threads: 8
    resources:
        mem_mb = 100000
    log: 
        expand("{logs}/{{individual}}/read_mapping.log", logs=config["log_dir"]),
    shell:
        """
        samtools fastq -@ {threads} -T MM,ML {input[0]} |\
        minimap2 -t {threads} -y -x map-hifi --MD --eqx --cs -Y -L -a -k 17 -K 10g {input[1]} - |\
        samtools view -@ {threads} -bh - |\
        samtools sort -@ {threads} - > {output}
        """

rule index_bam:
    input:
        expand("{bam_dir}/{{individual}}{extension}.bam", bam_dir = config['bam_dir'], extension = config['final_bam_extension'])
    output:
        expand("{bam_dir}/{{individual}}{extension}.bam.bai", bam_dir = config['bam_dir'], extension = config['final_bam_extension'])
    threads: 8
    log: expand("{logs}/{{individual}}/index.log", logs=config["log_dir"])
    shell:
        "samtools index -@ {threads} {input} {output} > {log} 2>&1"

