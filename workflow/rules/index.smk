rule index_vcf:
    input:
        expand("{vcf_dir}/{{vcf}}.vcf.gz", vcf_dir=config["vcf_dir"])
    output:
        expand("{vcf_dir}/{{vcf}}.vcf.gz.csi", vcf_dir=config["vcf_dir"])
    log:
        expand("{logs}/{{vcf}}.index.log", logs=config["log_dir"])
    shell:
        "bcftools index --threads {threads} {input} > {log} 2>&1"