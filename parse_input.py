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

    print(genotypes)
    return genotypes

def full_algorithm(file_path):
    f = open(file_path, 'r')
    lines = f.readlines()
    window_size = 2
    num_lines = len(lines)
    num_gens = len(lines[0])

    final_haps = []
    gen_len = 0
    full_gen_list = []
    full_final_haps = list()
    
    num_windows = (int)(num_lines/window_size)
    for i in range(0, num_windows): # break input into windows 
        for j in range(0, window_size):
            lines[window_size*i + j] = lines[window_size*i + j].rstrip().split()
        
        genotype_list = list()

        num_cols = len(lines[0])
        for c in range(0, num_cols):
            genotype = ""
            for r in range(0, window_size):
                genotype += lines[i*window_size+r][c]

            genotype_list.append(genotype)

        #print(genotype_list)
        gen_len = len(genotype_list[0])
        if i == 0:
            full_gen_list = genotype_list
        else:
            for x in range(0, len(genotype_list)):
                full_gen_list[x] += genotype_list[x]


        genotypes, haplotypes = EM(genotype_list)
        #print(genotypes)
        #print(haplotypes)

        for gen, phase_dict in genotypes.items():
            # determine the phase with the max probability
            maxProb = -1.0
            maxPhase = ('', '')
            for phase, prob in phase_dict.items():
                if prob > maxProb:
                    maxPhase = phase
                    maxProb = prob
            final_haps.append(maxPhase[0])
            final_haps.append(maxPhase[1])


###########################################################################################

    for j in range(0, num_lines - num_windows*window_size): #finish the final window
            lines[window_size*num_windows + j] = lines[window_size*num_windows + j].rstrip().split()

    genotype_list = list()

    for c in range(0, num_cols):
        genotype = ""
        for r in range(0, num_lines - num_windows*window_size):
            genotype += lines[num_windows*window_size+r][c]

        genotype_list.append(genotype)

    for x in range(0, len(genotype_list)):
        full_gen_list[x] += genotype_list[x]

    print(genotype_list)
    #print(full_gen_list)

    genotypes, haplotypes = EM(genotype_list)
    #print(genotypes)
    #print(haplotypes)
        
    for gen, phase_dict in genotypes.items():
        maxProb = -1.0
        maxPhase = ('', '')
        for phase, prob in phase_dict.items():
            if prob > maxProb:
                maxPhase = phase
                maxProb = prob
        final_haps.append(maxPhase[0])
        final_haps.append(maxPhase[1])


    print(final_haps)
    print(full_final_haps)
###########################################################################################

    # for i in range(0, gen_len):
    #     row = ''
    #     for h in final_haps:
    #         row += h[i] + ' '
    #     print(row) # replace with adding to file


def full2(file_path):
    f = open(file_path, 'r')
    input = f.readlines()
    window_size = 2
    num_snps = len(input)
    num_windows = (int)(num_snps/window_size)

    total_max_haps = []
    for i in range(0, num_windows):
        for r in range(0, window_size):
            input[window_size*i + r] = input[window_size*i + r].rstrip().split()
        genotype_list = []
        for c in range(0, len(input[0])):
            genotype = ''
            for r in range(0, window_size):
                genotype += input[i*window_size + r][c]
            genotype_list.append(genotype)

        print(genotype_list)
        genotypes, haplotypes = EM(genotype_list)
        max_haps = []
        for gen, phase_dict in genotypes.items():
            maxProb = -1.0
            maxPhase = ('', '')
            for phase, prob in phase_dict.items():
                if prob > maxProb:
                    maxPhase = phase
                    maxProb = prob
            if maxProb > 0:
                max_haps.extend(maxPhase)
        total_max_haps.append(max_haps)

    print(total_max_haps)
    

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
    print(genotypes)
    print(haplotypes)
    num_gens = len(genotypes)

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
            haplotypes[hap] = (probs[transient]/(2*num_gens), 0) 

    return (genotypes, haplotypes)


def main():
    #full2(sys.argv[1])
    #get_genotypes(sys.argv[1])

    geno, haplo = EM(['00', '10', '10'])
    print(geno)

    # genotype_list = get_genotypes(sys.argv[1])
    # gen_len = len(genotype_list[0])

    # genotypes, haplotypes = EM(genotype_list)
    
    # final_phases = []
    # final_haps = []
    # for gen, phase_dict in genotypes.items():
    #     maxProb = -1.0
    #     maxPhase = ('', '')
    #     for phase, prob in phase_dict.items():
    #         if prob > maxProb:
    #             maxPhase = phase
    #             maxProb = prob
    #     final_phases.append(maxPhase)
    #     final_haps.append(maxPhase[0])
    #     final_haps.append(maxPhase[1])

    # print(final_haps)


    # for i in range(0, gen_len):
    #     row = ''
    #     for h in final_haps:
    #         row += h[i] + ' '
    #     print(row) # replace with adding to file
    

if __name__ == "__main__":
    main()
