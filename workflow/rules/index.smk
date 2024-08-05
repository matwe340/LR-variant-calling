rule index_vcf:
    input:
        "{vcf}.vcf.gz"
    output:
        "{vcf}.vcf.gz.csi"
    log:
        expand("{logs}/{{vcf}}.index.log", logs=config["log_dir"])
    shell:
        "bcftools index --threads {threads} {input} > {log} 2>&1"