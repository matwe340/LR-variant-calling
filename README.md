# Variant-calling

Snakemake pipeline to do variant calling, that is, get from fastq files from the sequencing company to a vcf file.

## Installation

1. If you haven't yet, install [mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html) or [micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html).

2. Clone this repository and cd into it:

    ```
    git clone git@gitlab.gwdg.de:georg.manthey/variant-calling.git
    cd variant-calling
    ```

3. Create a new environment from the environment specs file:

    Using mamba:
    ```
    mamba env create -f environment.yml
    ```

    Using micromamba:
    ```
    micromamba env create -f environment.yml
    ```

## Usage

1. Copy the `config.yml.template` file to `config.yml` 

2. Adjust the paths to the genome and the raw fastq directory in the `config.yml` file.

    Depending on your setup and the sequencing technology used, change the path to the adapter sequences (the one provided uses illumina adapters provided by bbduk if you installed using conda, mamba or micromamba). Also adjust the optical duplicates distance depending on the sequencer used (e.g. 2500 for NovaSeq, 100 for MiSeq).

3. Create a chromosomes file from your reference genome:

    ```
    samtools faidx <reference.fasta> | cut -f 1 > chromosomes.txt
    ```

    And adjust the path in the `config.yml` (or place it in the data folder)

4. Create a individuals.txt file from your list of fastq files/sample sheet and adjust the read_name_pattern in the `config.yml` to fit your file names. The individuals.txt file needs to be a tab seperated file with 2 columns, the first one being the individual ids that should be in the final vcf and the second one should be the {id} part in the `read_name_pattern` of the `config.yml`. 

    For example:
    
    We have 2 individuals I1 and I2, with I1 being run on 2 lanes, and therefor resulting in files:
    
     - I1_L1_R1_001.fastq.gz
     - I1_L1_R2_001.fastq.gz
     - I1_L2_R1_001.fastq.gz
     - I1_L2_R2_001.fastq.gz
     - I2_L1_R1_001.fastq.gz
     - I2_L1_R2_001.fastq.gz

    Then the read_name_pattern should be:

    ```
    read_name_pattern: "{id}_R{read}_001.fastq.gz"
    ```

    And the individuals.txt file needs to look like:

    ```
    I1	I1_L1
    I1	I1_L2
    I2	I2_L1
    ```

    Then adjust the path in the `config.yml` file for the individuals file (or place it in the data folder).


    This allows for (almost) arbitrary files to be associated with any individual and a mapping to the actual individual ids early in the process. The trimmed fastq files will already contain combined sequences for an individual and are already named after the individual instead of after the files.

    This will likely be adjusted in the future though to allow for different read pair patterns for different files.

5. Run the pipeline.
    
    I suggest to open a screen or tmux window, as snakemake needs to run throughout the whole calculation in the background, but doesn't use many resources. 

    Activate the environment:

    with mamba:
    ```
    mamba activate variant-calling
    ```

    with micromamba:
    ```
    micromamba activate variant-calling
    ```

    Run
    ```
    snakemake --profile profile/default/
    ```

    in the variant-calling folder to start this pipeline on the rosa hpc cluster at the UOL. The profile sets sensible default resources for most use cases and limits the number of currently running jobs to 100, as there is a limit on number of concurrent jobs for a single group.

    If you want to only create bams and don't run variant calling, run
    ```
    snakemake --profile profile/default/ bams
    ```