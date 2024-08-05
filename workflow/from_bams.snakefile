import re

configfile: "config.yml"
localrules: filtered, index_bam, index_reference, merge_trimmed, retain_list

rule filtered:
    input:
        expand('{vcf_dir}/genome.IF-GF-MM2-RM.vcf.gz', vcf_dir = config["vcf_dir"]),
        expand('{vcf_dir}/genome.IF-GF-MM2-RM.vcf.gz.csi', vcf_dir = config["vcf_dir"])

include: "rules/scripts.smk"

include: "rules/call.smk"

include: "rules/filter.smk"