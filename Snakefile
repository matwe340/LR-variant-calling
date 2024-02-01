configfile: "config.yml"
localrules: filtered, unfiltered, index_bam, index_reference, merge_trimmed, index_vcf

def chromosome_vcfs(wildcards):
    with open(config["chromosome_file"], 'r') as f:
        chromosomes = [line.strip() for line in f]
    return expand("{vcf_dir}/{chromosome}.merged.raw.vcf.gz", vcf_dir = config["vcf_dir"], chromosome = chromosomes)

def chromosome_filtered_vcfs(wildcards):
    with open(config["chromosome_file"], 'r') as f:
        chromosomes = [line.strip() for line in f]
    return expand("{vcf_dir}/{chromosome}.merged.IF-GF-MM2.vcf.gz", vcf_dir = config["vcf_dir"], chromosome = chromosomes)


def get_individuals():  
    individuals = dict()
    with open(config["individual_file"], 'r') as f:
        for line in f:
            line = line.strip().split()
            if line[0] not in individuals:
                individuals[line[0]] = []
            individuals[line[0]].append(line[1])

    return individuals

rule filtered:
    input:
        chromosome_filtered_vcfs

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

rule index_vcf:
    input:
        "{vcf}.vcf.gz"
    output:
        "{vcf}.vcf.gz.csi"
    log:
        expand("{logs}/{{vcf}}.index.log", logs=config["log_dir"])
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

rule filter_qual_depth_missing_rpbz:
    input:
        expand("{vcf_dir}/{{chromosome}}.merged.raw.vcf.gz", vcf_dir = config["vcf_dir"])
    output:
        expand("{vcf_dir}/{{chromosome}}.merged.IF.vcf.gz", vcf_dir = config["vcf_dir"])
    log: expand("{logs}/{{chromosome}}/filter_qual_depth_missing.log", logs=config["log_dir"])
    run:
        sampn = int(shell("bcftools query -l {input} | wc -l", read=True))
        avgdp = int(shell("bcftools query -f '%DP\n' {input} | datamash median 1 | datamash round 1", read=True))
        dphi = 2 * avgdp

        shell(f"""bcftools view --types snps --threads {{threads}} -e "QUAL < 20 || INFO/DP > {dphi} || INFO/DP < {sampn} || MQ < 30 || RPBZ < -3 || RPBZ > 3" -Oz -o {output} {input} > {log} 2>&1""")

rule filter_gf:
    input:
        expand("{vcf_dir}/{{chromosome}}.merged.IF.vcf.gz", vcf_dir = config["vcf_dir"]),
        expand("{vcf_dir}/{{chromosome}}.merged.IF.vcf.gz.csi", vcf_dir = config["vcf_dir"])
    output:
        expand("{vcf_dir}/{{chromosome}}.merged.IF-GF.vcf.gz", vcf_dir = config["vcf_dir"])
    log: expand("{logs}/{{chromosome}}/filter_lcs.log", logs=config["log_dir"])
    shell:
        """bcftools +setGT -Oz -o {output} {input[0]} -- -t q -i "FMT/DP < {config[min_depth]}" -n "./." > {log} 2>&1"""

rule retain_list:
    input:
        expand("{vcf_dir}/{{chromosome}}.merged.IF-GF.vcf.gz", vcf_dir = config["vcf_dir"]),
        expand("{vcf_dir}/{{chromosome}}.merged.IF-GF.vcf.gz.csi", vcf_dir = config["vcf_dir"])
    output:
        expand("{vcf_dir}/{{chromosome}}_sample.stats", vcf_dir = config["vcf_dir"])
    shell:
        """echo -e "CHROM\tID\tnREF\tnALT\tnHET\tnTs\tnTv\tavgDP\tSingletons\tMissing_Sites" > {output}
        bcftools stats --threads {threads} -S- {input[0]} | grep 'PSC' | tr ' ' '_' | awk -v c={wildcards.chromosome} '{{OFS="\t"}}{{print c, $3,$4,$5,$6,$7,$8,$10,$11,$14}}' | sed '1,2d' >> {output}
        """

rule filter_genotype_missing_ind:
    input:
        expand("{vcf_dir}/{{chromosome}}.merged.IF-GF.vcf.gz", vcf_dir = config["vcf_dir"]),
        expand("{vcf_dir}/{{chromosome}}.merged.IF-GF.vcf.gz.csi", vcf_dir = config["vcf_dir"]),
        expand("{vcf_dir}/{{chromosome}}_sample.stats", vcf_dir = config["vcf_dir"])
    output:
        expand("{vcf_dir}/{{chromosome}}.merged.IF-GF-MM2.vcf.gz", vcf_dir = config["vcf_dir"])
    log: expand("{logs}/{{chromosome}}/filter_genotype_missing.log", logs=config["log_dir"])
    shell:
        """bcftools view --threads {threads} --samples-file {input[2]} --force-samples --Ou {input[0]} | bcftools view --min-ac 1 --threads {threads} -i 'F_MISSING<0.2' -Oz -o {output} > {log} 2>&1""" 

