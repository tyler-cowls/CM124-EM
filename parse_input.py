#from create_phases import phase
import create_phases
import sys

steady = 0
transient = 1

# obtain a list of genotypes from the input text file
# genotypes come out as a string of numbers with no spaces
def get_genotypes(file_path):
    f = open(file_path, 'r')
    snps = f.readlines()
    #print(snps)
    for i in range(0, len(snps)):
        snps[i] = snps[i].rstrip()
        snps[i] = snps[i].split()

    #print('Post-processing:')
    #print(snps)

    genotypes = list()

    for i in range(0, len(snps[0])):
        genotype = ""

        for j in range(0, len(snps)):
            genotype += snps[j][i]

        genotypes.append(genotype)

    return genotypes

# def full_algorithm(file_path):
#     f = open(file_path, 'r')
#     lines = f.readlines()
#     window_size = 50
#     num_lines = len(lines)
#     num_windows = 1+(int)(num_lines/window_size)
#     for i in range(0, num_windows): # break input into windows 
#         for j in range(0, window_size):
#             snps[window_size*i + j] = snps[window_size*i + j].rstrip().split()
        
#         genotypes = list()

#         l1 = len(snps[0])
#         for j in range(0, l1):
#             genotype = ""
#             l2 = len(snps)
#             for k in range(0, )


# build dictionary to map all haplotype pairs to their corresponding initial probabilities
def build_haplotype_pair_dictionary(haplotype_pair_list):
    num_pairs = len(haplotype_pair_list)
    haplotype_dict = {}
    initial_prob = 0.0#1.0 / num_pairs

    for i in range(0, num_pairs):
        haplotype_dict[haplotype_pair_list[i]] = initial_prob

    return haplotype_dict

# build dictionary that maps all genotypes to their corresponding haplotype dictionaries
# builds a unique haplotype dictionary
def build_dictionaries(genotypes):
    num_gens = len(genotypes)
    gen_dict = {}
    unique_hap_dict = {}
    unique_haps = set()

    for i in range(0, num_gens):
        (phases, haps) = create_phases.phase(genotypes[i])
        hap_dict = build_haplotype_pair_dictionary(phases)
        gen_dict[genotypes[i]] = hap_dict
        for hap in haps:
            unique_haps.add(hap)

    unique_haps = list(unique_haps)
    num_unique_haps = len(unique_haps)
    initial_hap_prob = 1.0 / num_unique_haps
    for i in range(0, num_unique_haps):
        unique_hap_dict[unique_haps[i]] = (initial_hap_prob, 0)

    return (gen_dict, unique_hap_dict)

# map all unique haplotypes to their corresponding initial probabilities
def build_unique_haplotype_dictionary(haplotype_dictionary):
    haps = set()
    for (key_gen, value_hap_pair) in haplotype_dictionary:
        hap_pair = value_hap_pair.key()

def EM(genotype_list):
    (genotypes, haplotypes) = build_dictionaries(genotype_list)

    for i in range(0, 10):
        for gen, phase_dict in genotypes.items():
            sumTempProb = 0
            for phase, prob in phase_dict.items():
                prob = haplotypes[phase[0]][steady]*haplotypes[phase[1]][steady]
                sumTempProb += prob
                genotypes[gen][phase] = prob

            for phase, prob in phase_dict.items():
                prob = prob/sumTempProb
                genotypes[gen][phase] = prob
                temp = haplotypes[phase[0]]
                haplotypes[phase[0]] = (temp[steady], temp[transient]+prob)
                temp = haplotypes[phase[1]]
                haplotypes[phase[1]] = (temp[steady], temp[transient]+prob) 

        # update haplotype probabilities
        for hap, probs in haplotypes.items():
            haplotypes[hap] = (probs[transient]/6, 0)

    return (genotypes, haplotypes)

def EM2(genotype_list):
    (genotypes, haplotypes) = build_dictionaries(genotype_list)

    # for each genotype -> (phase, prob)
    for gen, phase_dict in genotypes.items():
        # for each phase, compute phase probability
        # phase probability = P(h1)*P(h2) / (P(h1)*P(h2) + P(h3)*P(h4) + ...)
        sumTempProb = 0
        for phase, prob in phase_dict.items():
            prob = haplotypes[phase[0]][steady]*haplotypes[phase[1]][steady]
            sumTempProb += prob
        print(sumTempProb)
        for phase, prob in phase_dict.items():
            prob = prob / sumTempProb
            haplotypes[phase[0]] = (haplotypes[phase[0]][steady], haplotypes[phase[0]][transient]+prob)
            haplotypes[phase[1]] = (haplotypes[phase[1]][steady], haplotypes[phase[1]][transient]+prob)
            print(phase + ' - prob = ' + str(prob))


def main():
    #full_algorithm(sys.argv[1])
    #get_genotypes(sys.argv[1])


    genotype_list = get_genotypes(sys.argv[1])
    gen_len = len(genotype_list[0])

    (genotypes, haplotypes) = build_dictionaries(genotype_list)

    genotypes, haplotypes = EM(genotype_list)
    
    final_phases = []
    final_haps = []
    for gen, phase_dict in genotypes.items():
        maxProb = -1.0
        maxPhase = ('', '')
        for phase, prob in phase_dict.items():
            if prob > maxProb:
                maxPhase = phase
                maxProb = prob
        final_phases.append(maxPhase)
        final_haps.append(maxPhase[0])
        final_haps.append(maxPhase[1])


    for i in range(0, gen_len):
        row = ''
        for h in final_haps:
            row += h[i] + ' '
        print(row) # replace with adding to file
    

if __name__ == "__main__":
    main()
