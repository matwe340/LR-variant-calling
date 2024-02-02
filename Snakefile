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
    log: expand("{logs}/index_reference.log", logs=config["log_dir"])
    shell:
        "bwa-mem2 index -p {params[0]} {input[0]} > {log} 2>&1"

rule trim_reads:
    input:
        lambda wildcards: expand("{raw_fastq_dir}/" + config["read_name_pattern"], raw_fastq_dir = config["raw_fastq_dir"], read = [1, 2], id = wildcards.id),
    output:
        temp(expand("{fastq_trimmed_dir}/{{id}}_R{read}.trimmed.fastq.gz", fastq_trimmed_dir = config["fastq_trimmed_dir"], read = [1, 2]))
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

    return sorted(expand("{fastq_trimmed_dir}/{id}_R{{read}}.trimmed.fastq.gz", fastq_trimmed_dir = config["fastq_trimmed_dir"], id = individual_reads))

rule merge_trimmed:
    input:
        individual_trimmed
    output:
        expand("{fastq_trimmed_dir}/{{individual}}_R{{read}}.trimmed.all.fastq.gz", fastq_trimmed_dir = config["fastq_trimmed_dir"])
    shell:
        "cat {input} > {output}"

rule align:
    input:
        "output/genome.pac",
        expand("{fastq_trimmed_dir}/{{individual}}_R{read}.trimmed.all.fastq.gz", fastq_trimmed_dir = config["fastq_trimmed_dir"], read=[1, 2])
    params:
        genome_idx = "output/genome",
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
        "bwa-mem2 mem -t {threads} {params.genome_idx} {input[1]} {input[2]} 2> {log[0]} | samtools fixmate -@ {threads} -m - - 2> {log[1]} | samtools sort -@ {threads} -m {params.memory} -o {output} > {log[2]} 2>&1"

rule markdup:
    input:
        expand("{bam_dir}/{{individual}}.sorted.bam", bam_dir = config['bam_dir'])
    output:
        expand("{bam_dir}/{{individual}}.rmdup.bam", bam_dir = config['bam_dir'])
    log: expand("{logs}/{{individual}}/markdup.log", logs=config["log_dir"])
    threads: 4
    shell:
        "samtools markdup -@ {threads} -d {config[optical_dup_dist]} -S -r {input} {output} > {log} 2>&1"

rule index_bam:
    input:
        expand("{bam_dir}/{{individual}}.rmdup.bam", bam_dir = config['bam_dir'])
    output:
        expand("{bam_dir}/{{individual}}.rmdup.bam.bai", bam_dir = config['bam_dir'])
    log: expand("{logs}/{{individual}}/index.log", logs=config["log_dir"])
    shell:
        "samtools index -@ {threads} {input} {output} > {log} 2>&1"

rule call:
    input:
        expand("{bam_dir}/{{individual}}.rmdup.bam", bam_dir = config['bam_dir']),
        config["genome"],
        expand("{bam_dir}/{{individual}}.rmdup.bam.bai", bam_dir = config['bam_dir'])
    output:
        expand("{vcf_dir}/{{individual}}/{{chromosome}}.raw.vcf.gz", vcf_dir = config["vcf_dir"])
    threads: 2
    log: 
        expand("{logs}/{{individual}}/{{chromosome}}/mpielup.log", logs=config["log_dir"]),
        expand("{logs}/{{individual}}/{{chromosome}}/call.log", logs=config["log_dir"])
    shell:
        "bcftools mpileup --threads {threads} -q 20 -Q 20 -C 50 -Ou -r {wildcards.chromosome} -f {input[1]} {input[0]} -a \"AD,ADF,ADR,DP,SP\" 2> {log[0]} | bcftools call --threads {threads} --ploidy 2 -m -Oz -o {output} > {log[1]} 2>&1"

rule index_raw_vcf:
    input:
        expand("{vcf_dir}/{{individual}}/{{chromosome}}.raw.vcf.gz", vcf_dir = config['vcf_dir'])
    output:
        expand("{vcf_dir}/{{individual}}/{{chromosome}}.raw.vcf.gz.csi", vcf_dir = config['vcf_dir'])
    log:
        expand("{logs}/{{individual}}/{{chromosome}}/index.log", logs=config["log_dir"])
    shell:
        "bcftools index --threads {threads} {input} > {log} 2>&1"

def individual_vcfs(wildcards):
    individuals = get_individuals()

    vcf_list = expand("{vcf_dir}/{individual}/{chromosome}.raw.vcf.gz", vcf_dir = config["vcf_dir"], individual = individuals.keys(), chromosome = wildcards.chromosome)
    return vcf_list

def individual_vcf_indices(wildcards):
    individuals = get_individuals()

    vcf_list = expand("{vcf_dir}/{individual}/{chromosome}.raw.vcf.gz.csi", vcf_dir = config["vcf_dir"], individual = individuals.keys(), chromosome = wildcards.chromosome)
    return vcf_list

rule merge_vcf:
    input:
        vcf=individual_vcfs,
        index=individual_vcf_indices
    output:
        expand("{vcf_dir}/{{chromosome}}.merged.raw.vcf.gz", vcf_dir = config["vcf_dir"])
    threads: 2
    log: expand("{logs}/{{chromosome}}/merge.log", logs=config["log_dir"])
    shell:
        "bcftools merge --threads {threads} -Oz -o {output} {input.vcf} > {log} 2>&1"

