# Variant-calling

Snakemake pipeline to do variant calling, that is, get from fastq files from the sequencing company to a vcf file.

## Installation

1. If you haven't yet, install [mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html) or [micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html). ([conda](https://www.anaconda.com/) works as well, but I would recommend mamba/micromamba as they are a lot faster when installing)

2. Clone this repository and cd into it:

    ```bash
    git clone https://github.com/gmanthey/variant-calling.git
    cd variant-calling
    ```

3. Create a new environment from the environment specs file:

    Using mamba:
    ```bash
    mamba env create -f environment.yml
    ```

    Using micromamba:
    ```bash
    micromamba env create -f environment.yml
    ```

## Usage

1. Copy the `config.yml.template` file to `config.yml` 

2. Adjust the paths to the genome and the raw fastq directory in the `config.yml` file.

    Depending on your setup and the sequencing technology used, change the path to the adapter sequences (the one provided uses illumina adapters provided by bbduk if you installed using conda, mamba or micromamba). Also adjust the optical duplicates distance depending on the sequencer used (e.g. 2500 for NovaSeq, 100 for MiSeq).

3. Create a chromosomes file from your reference genome:

    ```bash
    samtools faidx <reference.fasta> 
    cut -f 1 <reference.fasta.fai> > chromosomes.txt
    ```
    If the reference contains many scaffolds, it may be more efficient to group scaffolds together (as this pipeline parallelizes by chromosomes). For this, change the chromosome file to contain a group name in the first column, followed by whitespace seperated sequence names. This format may also be combined with single sequence names for large sequences, i.e.:
    ```
    group1 seq1 seq2 seq3
    group2 seq4 seq5
    seq6
    seq7
    ```

    And adjust the path in the `config.yml` (or place it in the data folder)

4. Create a individuals.txt file from your list of fastq files/sample sheet. The individuals.txt file needs to be a tab seperated file with 2 columns, the first one being the individual ids that should be in the final vcf and the second one the path to the raw fastq file. Each line should only contain a single fastq file, but an individual can appear in multiple lines. The program automatically figures out R1 and R2 reads, as long as the string R1 or R2 appears in the file names.

    For example:
    
    We have 2 individuals I1 and I2, with I1 being run on 2 lanes, and therefor resulting in files:
    
     - I1_L1_R1_001.fastq.gz
     - I1_L1_R2_001.fastq.gz
     - I1_L2_R1_001.fastq.gz
     - I1_L2_R2_001.fastq.gz
     - I2_L1_R1_001.fastq.gz
     - I2_L1_R2_001.fastq.gz

    Then the individuals.txt file needs to look like:

    ```
    I1  I1_L1_R1_001.fastq.gz
    I1  I1_L1_R2_001.fastq.gz
    I1  I1_L2_R1_001.fastq.gz
    I1  I1_L2_R2_001.fastq.gz
    I2  I2_L1_R1_001.fastq.gz
    I2  I2_L1_R2_001.fastq.gz
    ```

    Then adjust the path in the `config.yml` file for the individuals file (or place it in the data folder).


    This allows for (almost) arbitrary files to be associated with any individual and a mapping to the actual individual ids early in the process. The trimmed fastq files will already contain combined sequences for an individual and are already named after the individual instead of after the files.

5. Run the pipeline.
    
    I suggest to open a screen or tmux window, as snakemake needs to run throughout the whole calculation in the background, but doesn't use many resources. 

    Activate the environment:

    with mamba:
    ```bash
    mamba activate variant-calling
    ```

    with micromamba:
    ```bash
    micromamba activate variant-calling
    ```

    Run
    ```bash
    snakemake --profile profile/default/
    ```

    in the variant-calling folder to start this pipeline on the rosa hpc cluster at the UOL. The profile sets sensible default resources for most use cases and limits the number of currently running jobs to 100, as there is a limit on number of concurrent jobs for a single group.

    If you want to only create bams and don't run variant calling, run
    ```bash
    snakemake --profile profile/default/ bams
    ```