configfile: "config.yml"
localrules: allsites_vcf

rule allsites_vcf:
    input:
        "allsites.vcf"

rule index:
    input: 
        config["genome"]
    output:
        "output/genome.idx"
    params:
        "output/genome"
    shell:
        "bwa-mem2 index -prefix {params[0]} {input[0]}"

rule align:
    input:
        "output/genome.idx",
        expand("{fastq_dir}/{{individual}}_R{read}.fastq.gz", fastq_dir = config["fastq_dir"], read = [1, 2])
    params:
        "output/genome"
    output:
        expand("{bam_dir}/{{individual}}.bam", bam_dir = config["bam_dir"])
    threads: 4
    log: "logs/{individual}/bwa.log"
    shell:
        "bwa-mem2 mem -t {threads} {params[0]} {input[1]} {input[2]} 2> {log[0]} | samtools view -b > {output}" 
    
