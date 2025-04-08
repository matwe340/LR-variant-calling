import os

include: "index.smk"

def raw_vcf_individual(wildcards):
    vcf_ro = f"{config['ro_ind_vcf_dir']}/{wildcards.individual}.raw.vcf.gz"
    if os.path.exists(vcf_ro):
        vcf_base = vcf_ro
    else:
        vcf_base = f"{config['vcf_dir']}/individuals/{wildcards.individual}.raw.vcf.gz"
    return [vcf_base, vcf_base + ".csi"]

rule filter_ind_quality:
    input:
        raw_vcf_individual
    output:
        temp(expand("{vcf_dir}/individuals/{{individual}}.QUAL.vcf.gz", vcf_dir = config["vcf_dir"]))
    log: expand("{logs}/{{individual}}/filter_qual.log", logs=config["log_dir"])
    run:
        shell(f"bcftools view --threads {{threads}} -e 'QUAL < 20' -Oz -o {{output}} {input[0]} > {{log}} 2>&1")

def individual_vcfs(wildcards):
    individuals = get_individuals()

    vcf_list = expand("{vcf_dir}/individuals/{individual}.QUAL.vcf.gz", vcf_dir = config["vcf_dir"], individual = individuals.keys())
    return vcf_list

def individual_vcf_indices(wildcards):
    individuals = get_individuals()

    vcf_list = expand("{vcf_dir}/individuals/{individual}.QUAL.vcf.gz.csi", vcf_dir = config["vcf_dir"], individual = individuals.keys())
    return vcf_list

def chromosome_filtered_vcfs(wildcards):
    chromosomes = get_chromosomes()
    return expand("{vcf_dir}/chromosomes/{chromosome}.IF-GF.vcf.gz", vcf_dir = config["vcf_dir"], chromosome = chromosomes.keys())

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

        shell(f"""bcftools view --types snps --threads {{threads}} -e "INFO/DP > {dphi} || INFO/DP < {sampn} || MQ < 30 || RPBZ < -3 || RPBZ > 3" -Oz -o {output} {input} > {log} 2>&1""")

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
        "results/retain.list"
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
        expand("results/retain.list", vcf_dir = config["vcf_dir"])
    output:
        expand("{vcf_dir}/genome.IF-GF-MM2.vcf.gz", vcf_dir = config["vcf_dir"])
    log: 
        expand("{logs}/filter_genotype_missing_samples.log", logs=config["log_dir"]),
        expand("{logs}/filter_genotype_missing_min.log", logs=config["log_dir"])
    shell:
        """bcftools view --threads {threads} --samples-file {input[2]} --force-samples -Ou {input[0]} 2> {log[0]} | bcftools view --min-ac 1 --threads {threads} -i 'F_MISSING<0.2' -Oz -o {output} > {log[1]} 2>&1""" 

rule sam_index_reference:
    input:
        config["genome"]
    output:
        "results/genome/genome.fai"
    log: expand("{logs}/sam_index_reference.log", logs=config["log_dir"])
    shell:
        "samtools faidx {input} --fai-idx {output} > {log} 2>&1"

rule filter_repeats:
    input:
        expand("{vcf_dir}/genome.IF-GF-MM2.vcf.gz", vcf_dir = config["vcf_dir"]),
        expand("{vcf_dir}/genome.IF-GF-MM2.vcf.gz.csi", vcf_dir = config["vcf_dir"]),
        "results/genome/genome.fai"
    output:
        expand("{vcf_dir}/genome.IF-GF-MM2-RM.vcf.gz", vcf_dir = config["vcf_dir"])
    log: expand("{logs}/filter_repeats.log", logs=config["log_dir"])
    shell:
        """bcftools view --threads {threads} -T <(bedtools complement -i {config[repeat_bed]} -g {input[2]}) -Oz -o {output} {input[0]} > {log} 2>&1"""

