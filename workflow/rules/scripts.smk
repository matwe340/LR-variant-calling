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