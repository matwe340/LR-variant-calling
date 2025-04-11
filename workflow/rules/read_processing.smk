##
import gzip
import os

def get_raw_fastq_files(wildcards):
    individuals = get_individuals()
    all_fastq_files = sum(individuals.values(), [])

    fastq_files = []
    for fastq_file in all_fastq_files:
        basename = os.path.basename(fastq_file)
        if basename in wildcards.run_id:
            fastq_files.append(fastq_file)

    return sorted(expand("{fastq_dir}/{fastq_file}", fastq_dir = config["raw_fastq_dir"], fastq_file = fastq_files), key=lambda x: x[::-1])

rule filter_and_trim_paired_reads:
    input:
        get_raw_fastq_files
    output:
        temp(expand("{fastq_trimmed_dir}/{{run_id}}_R{read}.trimmed.fastq.gz", fastq_trimmed_dir = config["fastq_trimmed_dir"], read=[1, 2])),
        "results/fastq_stats/{run_id}.json"
    log: expand("{logs}/{{run_id}}/trim.log", logs=config["log_dir"])
    threads: 4
    resources:
        mem_mb = 40000
    params:
        memory = "40G"
    shell:
        """fastp \
        --thread {threads} \
        --in1 {input[0]} --in2 {input[1]} \
        --out1 {output[0]} --out2 {output[1]} \
        --correction \
        --json {output[2]} \
        --html /dev/null
        """

def fastq_files_trimmed(wildcards):
    individuals = get_individuals()

    fastq_files_per_individual = individuals[wildcards.individual]

    run_ids = {}

    for fastq_file in fastq_files_per_individual:
        if fastq_file.endswith('gz'):
            file = gzip.open(os.path.join(config['raw_fastq_dir'], fastq_file), 'rt')
        else:
            file = open(os.path.join(config['raw_fastq_dir'], fastq_file), 'r')
        first_read = file.readline().strip().split()[0]
        file.close()
        first_read = first_read[1:].split('/')[0] # BGI names (and possibly others) have a /1 or /2 at the end
        if first_read not in run_ids:
            run_ids[first_read] = []
        run_ids[first_read].append(os.path.basename(fastq_file))

    for read_id, run_fastq_files in run_ids.items():
        if len(run_fastq_files) != 2:
            raise ValueError(f"Read id {read_id} does not have 2 files associate with it. (Found {', '.join(run_fastq_files)})")
    
    temp_file_ids = ['_'.join(sorted(run_fastq_files)) for run_fastq_files in run_ids.values()]

    return sorted(expand("{fastq_trimmed_dir}/{run_id}_R{read}.trimmed.fastq.gz", fastq_trimmed_dir = config["fastq_trimmed_dir"], run_id = temp_file_ids, read = wildcards.read))

rule merge_trimmed:
    input:
        fastq_files_trimmed
    output:
        expand("{fastq_trimmed_dir}/{{individual}}_R{{read}}.trimmed.all.fastq.gz", fastq_trimmed_dir = config["fastq_trimmed_dir"])
    shell:
        "cat {input} > {output}"

def get_summary_files(wildcards):
    individuals = get_individuals()

    temp_file_ids = []
    for individual in individuals:
        fastq_files_per_individual = individuals[individual]

        run_ids = {}

        for fastq_file in fastq_files_per_individual:
            if fastq_file.endswith('gz'):
                file = gzip.open(os.path.join(config['raw_fastq_dir'], fastq_file), 'rt')
            else:
                file = open(os.path.join(config['raw_fastq_dir'], fastq_file), 'r')
            first_read = file.readline().strip().split()[0]
            file.close()
            first_read = first_read[1:].split('/')[0] # BGI names (and possibly others) have a /1 or /2 at the end
            if first_read not in run_ids:
                run_ids[first_read] = []
            run_ids[first_read].append(os.path.basename(fastq_file))

        for read_id, run_fastq_files in run_ids.items():
            if len(run_fastq_files) != 2:
                raise ValueError(f"Read id {read_id} does not have 2 files associate with it. (Found {', '.join(run_fastq_files)})")
    
        temp_file_ids.extend(['_'.join(sorted(run_fastq_files)) for run_fastq_files in run_ids.values()])

    return expand("results/fastq_stats/{run_id}.summary", run_id = temp_file_ids)

rule summarize_fastq_stats:
    input:
        get_summary_files
    output:
        "results/fastq_stats/runs_fastq_summary"
    shell:
        "cat {input} > {output}"

rule parse_json:
    input:
        "results/fastq_stats/{run_id}.json"
    output:
        "results/fastq_stats/{run_id}.summary"
    shell:
        """
        python helpers/parse_json.py \
        --fastp-json {input} \
        --o results/fastq_stats \
        -s {wildcards.run_id}
        """




