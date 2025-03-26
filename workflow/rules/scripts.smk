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
            if len(line) == 1:
                individuals[line[0]].append(line[0])
            else:
                individuals[line[0]].append(line[1])

    return individuals

def get_bam_file(individual):
    if os.path.exists(f"{config["ro_bam_dir"]}/{individual}{config['final_bam_extension']}.bam"):
        bam_prefix = f"{config["ro_bam_dir"]}/{individual}{config['final_bam_extension']}"
    else:
        bam_prefix = f"{config['bam_dir']}/{individual}{config['final_bam_extension']}"

    bam_file = bam_prefix + ".bam"

    if os.path.exists(bam_prefix + ".bam") and os.path.exists(bam_prefix + ".bai"):
            return [bam_file, bam_prefix + ".bai"]

    return [bam_file, bam_prefix + ".bam.bai"]