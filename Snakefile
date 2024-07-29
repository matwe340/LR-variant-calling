import re

configfile: "config.yml"
localrules: filtered, index_bam, index_reference, merge_trimmed, index_raw_vcf, retain_list

def get_chromosomes():
    chromosomes = {}
    with open(config["chromosome_file"], 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            if len(line) == 1:
                chromosomes[line[0]] = line
            else:
                chromosomes[line[0]] = line[1:]

    return chromosomes

def chromosome_filtered_vcfs(wildcards):
    chromosomes = get_chromosomes()
    return expand("{vcf_dir}/chromosomes/{chromosome}.IF-GF.vcf.gz", vcf_dir = config["vcf_dir"], chromosome = chromosomes.keys())

def get_individuals():  
    individuals = dict()
    with open(config["individual_file"], 'r') as f:
        for line in f:
            line = line.strip()
            if line == '':
                continue
            line = line.split()
            if line[0] not in individuals:
                individuals[line[0]] = []
            individuals[line[0]].append(line[1])

    return individuals

rule filtered:
    input:
        expand('{vcf_dir}/genome.IF-GF-MM2-RM.vcf.gz', vcf_dir = config["vcf_dir"]),
        expand('{vcf_dir}/genome.IF-GF-MM2-RM.vcf.gz.csi', vcf_dir = config["vcf_dir"])

def individual_bams(wildcards):
    individuals = get_individuals()
    bam_list = expand("{bam_dir}/{individual}{extension}.bam", bam_dir = config["bam_dir"], individual = individuals.keys(), extension=config["final_bam_extension"])
    bam_index_list = expand("{bam_dir}/{individual}{extension}.bai", bam_dir = config["bam_dir"], individual = individuals.keys(), extension=config["final_bam_extension"])
    return bam_list + bam_index_list

rule bams:
    input:
        individual_bams

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

rule sam_index_reference:
    input:
        config["genome"]
    output:
        "output/genome.fai"
    log: expand("{logs}/sam_index_reference.log", logs=config["log_dir"])
    shell:
        "samtools faidx {input} --fai-idx {output} > {log} 2>&1"

def get_reads(wildcards):
    individual_reads = get_individuals()
    all_reads = sum(individual_reads.values(), start=[])
    reads = sorted([read for read in all_reads if read.startswith(wildcards.id)])

    return expand("{fastq_dir}/{read}", fastq_dir = config["raw_fastq_dir"], read = reads)

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
        read_id_pos = re.search(r"R[12]", read).start()
        individual_id = read[:read_id_pos]
        if individual_id not in individual_ids:
            individual_ids[individual_id] = []
        individual_ids[individual_id].append(read)
    
    for id, reads in individual_ids.items():
        if len(reads) > 2:
            raise ValueError(f"Individual {id} has more than 2 reads ({reads})")

    return sorted(expand("{fastq_trimmed_dir}/{id}_R{{read}}.trimmed.fastq.gz", fastq_trimmed_dir = config["fastq_trimmed_dir"], id = individual_ids.keys()))

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
        expand("{bam_dir}/{{individual}}{extension}.bam", bam_dir = config['bam_dir'], extension = config['final_bam_extension'])
    log: expand("{logs}/{{individual}}/markdup.log", logs=config["log_dir"])
    threads: 4
    shell:
        "samtools markdup -@ {threads} -d {config[optical_dup_dist]} -S -r {input} {output} > {log} 2>&1"

rule index_bam:
    input:
        expand("{bam_dir}/{{individual}}{extension}.bam", bam_dir = config['bam_dir'], extension = config['final_bam_extension'])
    output:
        expand("{bam_dir}/{{individual}}{extension}.bai", bam_dir = config['bam_dir'], extension = config['final_bam_extension'])
    log: expand("{logs}/{{individual}}/index.log", logs=config["log_dir"])
    shell:
        "samtools index -@ {threads} {input} {output} > {log} 2>&1"

rule individual_file:
    input:
        config["individual_file"]
    output:
        temp("output/{individual}.txt")
    run:
        with open(output[0], "w") as f:
            f.write(wildcards.individual)

def bam_index_file(wildcards):
    if os.path.exists(f"{config['bam_dir']}/{wildcards.individual}{config['final_bam_extension']}.bai"):
        return f"{config['bam_dir']}/{wildcards.individual}{config['final_bam_extension']}.bai"
    else:
        return f"{config['bam_dir']}/{wildcards.individual}{config['final_bam_extension']}.bam.bai"

rule call:
    input:
        expand("{bam_dir}/{{individual}}{extension}.bam", bam_dir = config['bam_dir'], extension = config['final_bam_extension']),
        config["genome"],
        bam_index_file,
        "output/{individual}.txt"
    output:
        temp(f"{config['vcf_dir']}/{{individual}}.raw.unnamed.vcf.gz")
    threads: 2
    log: 
        expand("{logs}/{{individual}}/mpielup.log", logs=config["log_dir"]),
        expand("{logs}/{{individual}}/call.log", logs=config["log_dir"])
    shell:
        "bcftools mpileup --threads {threads} -q 20 -Q 20 -Ou -s {wildcards.individual} --ignore-RG -f {input[1]} {input[0]} -a \"AD,ADF,ADR,DP,SP\" 2> {log[0]} | bcftools call --threads {threads} --ploidy 2 -m -Oz -o {output} > {log[1]} 2>&1"

rule rename_individual:
    input:
        expand("{vcf_dir}/{{individual}}.raw.unnamed.vcf.gz", vcf_dir=config["vcf_dir"]),
        "output/{individual}.txt"
    output:
        expand("{vcf_dir}/{{individual}}.raw.vcf.gz", vcf_dir=config["vcf_dir"])
    log:
        expand("{logs}/{{individual}}/rename.log", logs=config["log_dir"])
    threads: 2
    shell:
        "bcftools reheader --threads {threads} -s {input[1]} {input[0]} -Oz -o {output} > {log} 2>&1"

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

    vcf_list = expand("{vcf_dir}/{individual}.raw.vcf.gz", vcf_dir = config["vcf_dir"], individual = individuals.keys())
    return vcf_list

def individual_vcf_indices(wildcards):
    individuals = get_individuals()

    vcf_list = expand("{vcf_dir}/{individual}.raw.vcf.gz.csi", vcf_dir = config["vcf_dir"], individual = individuals.keys())
    return vcf_list

rule merge_vcf_into_chromosomes:
    input:
        vcf=individual_vcfs,
        index=individual_vcf_indices
    output:
        expand("{vcf_dir}/chromosomes/{{chromosome}}.allsites.vcf.gz", vcf_dir = config["vcf_dir"])
    threads: 2
    log: expand("{logs}/{{chromosome}}/merge.log", logs=config["log_dir"])
    run:
        chromosomes = ','.join(get_chromosomes()[wildcards.chromosome])
        shell(f"bcftools merge --threads {{threads}} -r {chromosomes} -Oz -o {{output}} {{input.vcf}} > {{log}} 2>&1")

rule filter_qual_depth_missing_rpbz:
    input:
        expand("{vcf_dir}/chromosomes/{{chromosome}}.allsites.vcf.gz", vcf_dir = config["vcf_dir"])
    output:
        expand("{vcf_dir}/chromosomes/{{chromosome}}.IF.vcf.gz", vcf_dir = config["vcf_dir"])
    log: expand("{logs}/{{chromosome}}/filter_qual_depth_missing.log", logs=config["log_dir"])
    run:
        sampn = int(shell("bcftools query -l {input} | wc -l", read=True))
        avgdp = int(shell("bcftools query -f '%DP\n' {input} | datamash median 1 | datamash round 1", read=True))
        dphi = 2 * avgdp

        shell(f"""bcftools view --types snps --threads {{threads}} -e "QUAL < 20 || INFO/DP > {dphi} || INFO/DP < {sampn} || MQ < 30 || RPBZ < -3 || RPBZ > 3" -Oz -o {output} {input} > {log} 2>&1""")

rule filter_gf:
    input:
        expand("{vcf_dir}/chromosomes/{{chromosome}}.IF.vcf.gz", vcf_dir = config["vcf_dir"]),
        expand("{vcf_dir}/chromosomes/{{chromosome}}.IF.vcf.gz.csi", vcf_dir = config["vcf_dir"])
    output:
        expand("{vcf_dir}/chromosomes/{{chromosome}}.IF-GF.vcf.gz", vcf_dir = config["vcf_dir"])
    log: expand("{logs}/{{chromosome}}/filter_lcs.log", logs=config["log_dir"])
    shell:
        """bcftools +setGT -Oz -o {output} {input[0]} -- -t q -i "FMT/DP < {config[min_depth]}" -n "./." > {log} 2>&1"""

rule merge_chromosome_vcf:
    input:
        chromosome_filtered_vcfs
    output:
        expand("{vcf}/genome.IF-GF.vcf.gz", vcf = config["vcf_dir"])
    log: expand("{logs}/merge_all.log", logs=config["log_dir"])
    shell:
        "bcftools concat -Oz -o {output} {input} > {log} 2>&1"

rule sample_stats:
    input:
        expand("{vcf_dir}/genome.IF-GF.vcf.gz", vcf_dir = config["vcf_dir"]),
        expand("{vcf_dir}/genome.IF-GF.vcf.gz.csi", vcf_dir = config["vcf_dir"])
    output:
        expand("{vcf_dir}/sample.stats", vcf_dir = config["vcf_dir"])
    shell:
        """echo -e "ID\tnREF\tnALT\tnHET\tnTs\tnTv\tavgDP\tSingletons\tMissing_Sites" > {output}
        bcftools stats --threads {threads} -S- {input[0]} | grep 'PSC' | tr ' ' '_' | awk '{{OFS="\t"}}{{print $3,$4,$5,$6,$7,$8,$10,$11,$14}}' | sed '1,2d' >> {output}
        """

rule retain_list:
    input:
        expand("{vcf_dir}/sample.stats", vcf_dir = config["vcf_dir"])
    output:
        "data/retain.list"
    run:
        with open(input[0], 'r') as f:
            f.readline()
            individuals = [line.strip().split()[0] for line in f]
        with open(output[0], 'w') as f:
            f.write('\n'.join(individuals))        
    

rule filter_genotype_missing_ind:
    input:
        expand("{vcf_dir}/genome.IF-GF.vcf.gz", vcf_dir = config["vcf_dir"]),
        expand("{vcf_dir}/genome.IF-GF.vcf.gz.csi", vcf_dir = config["vcf_dir"]),
        expand("data/retain.list", vcf_dir = config["vcf_dir"])
    output:
        expand("{vcf_dir}/genome.IF-GF-MM2.vcf.gz", vcf_dir = config["vcf_dir"])
    log: 
        expand("{logs}/filter_genotype_missing_samples.log", logs=config["log_dir"]),
        expand("{logs}/filter_genotype_missing_min.log", logs=config["log_dir"])
    shell:
        """bcftools view --threads {threads} --samples-file {input[2]} --force-samples -Ou {input[0]} 2> {log[0]} | bcftools view --min-ac 1 --threads {threads} -i 'F_MISSING<0.2' -Oz -o {output} > {log[1]} 2>&1""" 

rule filter_repeats:
    input:
        expand("{vcf_dir}/genome.IF-GF-MM2.vcf.gz", vcf_dir = config["vcf_dir"]),
        expand("{vcf_dir}/genome.IF-GF-MM2.vcf.gz.csi", vcf_dir = config["vcf_dir"]),
        "output/genome.fai"
    output:
        expand("{vcf_dir}/genome.IF-GF-MM2-RM.vcf.gz", vcf_dir = config["vcf_dir"])
    log: expand("{logs}/filter_repeats.log", logs=config["log_dir"])
    shell:
        """bcftools view --threads {threads} -T <(bedtools complement -i {config[repeat_bed]} -g {input[2]}) -Oz -o {output} {input[0]} > {log} 2>&1"""

