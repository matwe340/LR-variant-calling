import gzip
import os

def individual_bams(wildcards):
    individuals = get_individuals()
    return flatten([get_bam_file(individual) for individual in individuals.keys()])

def get_raw_fastq_files(wildcards):
    individuals = get_individuals()
    all_fastq_files = sum(individuals.values(), [])

    fastq_files = []
    for fastq_file in all_fastq_files:
        basename = os.path.basename(fastq_file)
        if basename in wildcards.run_id:
            fastq_files.append(fastq_file)

    return sorted(expand("{fastq_dir}/{fastq_file}", fastq_dir = config["raw_fastq_dir"], fastq_file = fastq_files), key=lambda x: x[::-1])

rule bams:
    input:
        individual_bams

rule index_reference:
    input: 
        config["genome"]
    output:
        "results/genome/genome.pac"
    params:
        "results/genome/genome"
    log: expand("{logs}/index_reference.log", logs=config["log_dir"])
    shell:
        "bwa-mem2 index -p {params[0]} {input[0]} > {log} 2>&1"

def trimmed_fastq_individual(wildcards):
    fastq_ro = expand("{ro_fastq_trimmed_dir}/{individual}_R{read}.trimmed.all.fastq.gz", ro_fastq_trimmed_dir = config["ro_fastq_trimmed_dir"], individual = wildcards.individual, read = [1, 2])
    if all([os.path.exists(fastq_file) for fastq_file in fastq_ro]):
        return fastq_ro
    else:
        return expand("{fastq_trimmed_dir}/{individual}_R{read}.trimmed.all.fastq.gz", fastq_trimmed_dir = config["fastq_trimmed_dir"], individual = wildcards.individual, read = [1, 2])

rule align:
    input:
        "results/genome/genome.pac",
        trimmed_fastq_individual
    params:
        genome_idx = "results/genome/genome",
        memory = "8G"
    output:
        temp(expand("{bam_dir}/{{individual}}.sorted.bam", bam_dir = config["bam_dir"]))
    threads: 8
    resources:
        mem_mb = 100000
    log: 
        expand("{logs}/{{individual}}/bwa.log", logs=config["log_dir"]),
        expand("{logs}/{{individual}}/fixmate.log", logs=config["log_dir"]),
        expand("{logs}/{{individual}}/sort.log", logs=config["log_dir"])
    shell:
        "bwa-mem2 mem -t {threads} {params.genome_idx} {input[1]} {input[2]} 2> {log[0]} | samtools fixmate -@ {threads} -u -m - - 2> {log[1]} | samtools sort -@ {threads} -m {params.memory} -o {output} > {log[2]} 2>&1"

rule markdup:
    input:
        expand("{bam_dir}/{{individual}}.sorted.bam", bam_dir = config['bam_dir'])
    output:
        expand("{bam_dir}/{{individual}}{extension}.bam", bam_dir = config['bam_dir'], extension = config['final_bam_extension'])
    log: expand("{logs}/{{individual}}/markdup.log", logs=config["log_dir"])
    threads: 4
    shell:
        "samtools markdup -@ {threads} -S -r {input} {output} > {log} 2>&1"

rule index_bam:
    input:
        expand("{bam_dir}/{{individual}}{extension}.bam", bam_dir = config['bam_dir'], extension = config['final_bam_extension'])
    output:
        expand("{bam_dir}/{{individual}}{extension}.bam.bai", bam_dir = config['bam_dir'], extension = config['final_bam_extension'])
    log: expand("{logs}/{{individual}}/index.log", logs=config["log_dir"])
    shell:
        "samtools index -@ {threads} {input} {output} > {log} 2>&1"