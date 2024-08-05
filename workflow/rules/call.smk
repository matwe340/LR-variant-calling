import os

def bam_index_file(wildcards):
    if os.path.exists(f"{config['bam_dir']}/{wildcards.individual}{config['final_bam_extension']}.bam"):
        if os.path.exists(f"{config['bam_dir']}/{wildcards.individual}{config['final_bam_extension']}.bai"):
            return f"{config['bam_dir']}/{wildcards.individual}{config['final_bam_extension']}.bai"
        else:
            return f"{config['bam_dir']}/{wildcards.individual}{config['final_bam_extension']}.bam.bai"
    else:
        return f"{config['bam_dir']}/{wildcards.individual}{config['final_bam_extension']}.bam.bai"

rule individual_file:
    input:
        config["individual_file"]
    output:
        temp("results/{individual}.txt")
    run:
        with open(output[0], "w") as f:
            f.write(wildcards.individual)

rule call:
    input:
        expand("{bam_dir}/{{individual}}{extension}.bam", bam_dir = config['bam_dir'], extension = config['final_bam_extension']),
        config["genome"],
        bam_index_file,
        "results/{individual}.txt"
    output:
        temp(f"{config['vcf_dir']}/{{individual}}.raw.unnamed.vcf.gz")
    threads: 2
    log: 
        expand("{logs}/{{individual}}/mpielup.log", logs=config["log_dir"]),
        expand("{logs}/{{individual}}/call.log", logs=config["log_dir"])
    shell:
        "bcftools mpileup --threads {threads} -q 20 -Q 20 -Ou -s {wildcards.individual} --ignore-RG -f {input[1]} {input[0]} -a \"AD,ADF,ADR,DP,SP\" 2> {log[0]} | bcftools call --threads {threads} -a \"GQ\" --ploidy 2 -m -Oz -o {output} > {log[1]} 2>&1"

rule rename_individual:
    input:
        expand("{vcf_dir}/{{individual}}.raw.unnamed.vcf.gz", vcf_dir=config["vcf_dir"]),
        "results/{individual}.txt"
    output:
        expand("{vcf_dir}/{{individual}}.raw.vcf.gz", vcf_dir=config["vcf_dir"])
    log:
        expand("{logs}/{{individual}}/rename.log", logs=config["log_dir"])
    threads: 2
    shell:
        "bcftools reheader --threads {threads} -s {input[1]} {input[0]} -o {output} > {log} 2>&1"

include: "index.smk"