import gzip
import os

def individual_bams(wildcards):
    individuals = get_individuals()
    bam_list = expand("{bam_dir}/{individual}{extension}.bam", bam_dir = config["bam_dir"], individual = individuals.keys(), extension=config["final_bam_extension"])
    bam_index_list = expand("{bam_dir}/{individual}{extension}.bai", bam_dir = config["bam_dir"], individual = individuals.keys(), extension=config["final_bam_extension"])
    return bam_list + bam_index_list

def get_reads(wildcards):
    individuals = get_individuals()
    all_reads = sum(individuals.values(), [])

    reads = []
    for read in all_reads:
        basename = os.path.basename(read)
        if basename in wildcards.id:
            reads.append(read)

    return sorted(expand("{fastq_dir}/{read}", fastq_dir = config["raw_fastq_dir"], read = reads), key=lambda x: x[::-1])

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

rule trim_paired_reads:
    input:
        get_reads
    output:
        temp(expand("{fastq_trimmed_dir}/{{id}}_R{read}.trimmed.fastq.gz", fastq_trimmed_dir = config["fastq_trimmed_dir"], read=[1, 2]))
    log: expand("{logs}/{{id}}/trim.log", logs=config["log_dir"])
    threads: 4
    resources:
        mem_mb = 40000
    params:
        memory = "40G"
    shell:
        "bbduk.sh t={threads} -Xmx{params.memory} overwrite=true in={input[0]} in2={input[1]} out={output[0]} out2={output[1]} ref={config[adapters]} ktrim=r k=23 mink=25 hdist=1 tpe tbo > {log[0]} 2>&1"

def individual_trimmed(wildcards):
    individuals = get_individuals()

    individual_reads = individuals[wildcards.individual]

    individual_ids = {}

    for read in individual_reads:
        if read.endswith('gz'):
            file = gzip.open(os.path.join(config['raw_fastq_dir'], read), 'rt')
        else:
            file = open(os.path.join(config['raw_fastq_dir'], read), 'r')
        read1 = file.readline().strip().split()[0]
        file.close()
        read1 = read1[1:].split('/')[0] # BGI names (and possibly others) have a /1 or /2 at the end
        if read1 not in individual_ids:
            individual_ids[read1] = []
        individual_ids[read1].append(os.path.basename(read))

    for read_id, read_files in individual_ids.items():
        if len(read_files) != 2:
            raise ValueError(f"Read id {read_id} does not have 2 files associate with it. (Found {', '.join(read_files)})")
    
    file_ids = ['_'.join(sorted(read_files)) for read_files in individual_ids.values()]

    return sorted(expand("{fastq_trimmed_dir}/{id}_R{read}.trimmed.fastq.gz", fastq_trimmed_dir = config["fastq_trimmed_dir"], id = file_ids, read = wildcards.read))

rule merge_trimmed:
    input:
        individual_trimmed
    output:
        expand("{fastq_trimmed_dir}/{{individual}}_R{{read}}.trimmed.all.fastq.gz", fastq_trimmed_dir = config["fastq_trimmed_dir"])
    shell:
        "cat {input} > {output}"

rule align:
    input:
        "results/genome/genome.pac",
        expand("{fastq_trimmed_dir}/{{individual}}_R{read}.trimmed.all.fastq.gz", fastq_trimmed_dir = config["fastq_trimmed_dir"], read=[1, 2])
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
