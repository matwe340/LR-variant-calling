configfile: "config.yml"
localrules: unfiltered, index_bam, index_reference, merge_trimmed, index_raw_vcf

def chromosome_vcfs(wildcards):
    with open(config["chromosome_file"], 'r') as f:
        chromosomes = [line.strip() for line in f]
    return expand("{vcf_dir}/{chromosome}.merged.raw.vcf.gz", vcf_dir = config["vcf_dir"], chromosome = chromosomes)

def get_individuals():  
    individuals = dict()
    with open(config["individual_file"], 'r') as f:
        for line in f:
            line = line.strip().split()
            if line[0] not in individuals:
                individuals[line[0]] = []
            individuals[line[0]].append(line[1])

    return individuals

rule unfiltered:
    input:
        chromosome_vcfs

rule index_reference:
    input: 
        config["genome"]
    output:
        "output/genome.pac"
    params:
        "output/genome"
    shell:
        "bwa-mem2 index -p {params[0]} {input[0]}"

rule trim_reads:
    input:
        lambda wildcards: expand("{raw_fastq_dir}/" + config["read_name_pattern"], raw_fastq_dir = config["raw_fastq_dir"], read = [1, 2], id = wildcards.id),
    output:
        temp(expand("{fastq_trimmed_dir}/{{id}}.trimmed.fastq.gz", fastq_trimmed_dir = config["fastq_trimmed_dir"]))
    log: "logs/{id}/trim.log"
    threads: 4
    resources:
        mem_mb = 40000
    params:
        memory = "40G"
    shell:
        "bbduk.sh t={threads} -Xmx{params.memory} overwrite=true in={input[0]} in2={input[1]} out={output[0]} ref={config[adapters]} ktrim=r k=23 mink=25 hdist=1 tpe tbo > {log[0]} 2>&1"

def individual_trimmed(wildcards):
    individuals = get_individuals()

    individual_reads = individuals[wildcards.individual]

    return expand("{fastq_trimmed_dir}/{id}.trimmed.fastq.gz", fastq_trimmed_dir = config["fastq_trimmed_dir"], id = individual_reads)

rule merge_trimmed:
    input:
        individual_trimmed
    output:
        expand("{fastq_trimmed_dir}/{{individual}}.trimmed.all.fastq.gz", fastq_trimmed_dir = config["fastq_trimmed_dir"])
    shell:
        "cat {input} > {output}"

rule align:
    input:
        "output/genome.pac",
        expand("{fastq_trimmed_dir}/{{individual}}.trimmed.all.fastq.gz", fastq_trimmed_dir = config["fastq_trimmed_dir"])
    params:
        genome_idx = "output/genome",
        memory = "8G"
    output:
        temp(expand("{bam_dir}/{{individual}}.sorted.bam", bam_dir = config["bam_dir"]))
    threads: 8
    resources:
        mem_mb = 100000
    log: "logs/{individual}/bwa.log"
    shell:
        "bwa-mem2 mem -t {threads} {params.genome_idx} {input[1]} 2> {log[0]} | samtools sort -@ {threads} -m {params.memory} -o {output}" 
    
rule markdup:
    input:
        expand("{bam_dir}/{{individual}}.sorted.bam", bam_dir = config['bam_dir'])
    output:
        expand("{bam_dir}/{{individual}}.rmdup.bam", bam_dir = config['bam_dir'])
    log: "logs/{individual}/markdup.log"
    threads: 4
    shell:
        "samtools markdup -@ {threads} -S -r {input} {output} 2> {log}"

rule index_bam:
    input:
        expand("{bam_dir}/{{individual}}.rmdup.bam", bam_dir = config['bam_dir'])
    output:
        expand("{bam_dir}/{{individual}}.rmdup.bam.bai", bam_dir = config['bam_dir'])
    shell:
        "samtools index -@ {threads} {input} {output}"

rule call:
    input:
        expand("{bam_dir}/{{individual}}.rmdup.bam", bam_dir = config['bam_dir']),
        config["genome"],
        expand("{bam_dir}/{{individual}}.rmdup.bam.bai", bam_dir = config['bam_dir'])
    output:
        expand("{vcf_dir}/{{individual}}/{{chromosome}}.raw.vcf.gz", vcf_dir = config["vcf_dir"])
    threads: 2
    shell:
        "bcftools mpileup --threads {threads} -q 20 -Q 20 -C 50 -Ou -r {wildcards.chromosome} -f {input[1]} {input[0]} -a \"AD,ADF,ADR,DP,SP\" | bcftools call --threads {threads} --ploidy 2 -m -Oz -o {output}"

rule index_raw_vcf:
    input:
        expand("{vcf_dir}/{{individual}}/{{chromosome}}.raw.vcf.gz", vcf_dir = config['vcf_dir'])
    output:
        expand("{vcf_dir}/{{individual}}/{{chromosome}}.raw.vcf.gz.csi", vcf_dir = config['vcf_dir'])
    shell:
        "bcftools index --threads {threads} {input}"

def individual_vcfs(wildcards):
    individuals = get_individuals()

    vcf_list = expand("{vcf_dir}/{individual}/{chromosome}.raw.vcf.gz", vcf_dir = config["vcf_dir"], individual = individuals.keys(), chromosome = wildcards.chromosome)
    vcf_list.extend(expand("{vcf_dir}/{individual}/{chromosome}.raw.vcf.gz.csi", vcf_dir = config["vcf_dir"], individual = individuals.keys(), chromosome = wildcards.chromosome))
    return vcf_list

rule merge_vcf:
    input:
        individual_vcfs
    output:
        expand("{vcf_dir}/{{chromosome}}.merged.raw.vcf.gz", vcf_dir = config["vcf_dir"])
    threads: 2
    shell:
        "bcftools merge --threads {threads} -Oz -o {output} {input[:len(inout)/2]}"
