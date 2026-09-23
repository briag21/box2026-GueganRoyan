import random

def simulate_reads_of_fixed_size(genome, k, seed=None):
    '''
    Generate a set of single-strain errorless reads of lenght k that cover the
    genome entirely
    '''
    if seed != None:
        random.seed(seed)

    i = 0
    reads = []
    while i < len(genome) - k:
        reads.append(genome[i:i+k])
        i += random.randint(1, k)
    reads.append(genome[-k:])
    return reads

