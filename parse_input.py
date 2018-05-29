from create_phases import phase

# obtain a list of genotypes from the input text file
# genotypes come out as a string of numbers with no spaces
def get_genotypes(file_path):
    f = open(file_path, 'r')
    snps = f.readlines()
    for i in range(0, len(snps)):
        snps[i] = snps[i].rstrip()
        snps[i] = snps[i].split()

    genotypes = list()

    for i in range(0, len(snps[0])):
        genotype = ""

        for j in range(0, len(snps)):
            genotype += snps[j][i]

        genotypes.append(genotype)

    return genotypes

# build dictionary to map all haplotype pairs to their corresponding initial probabilities
def build_haplotype_pair_dictionary(haplotype_pair_list):
    num_pairs = len(haplotype_pair_list)
    haplotype_dict = {}
    initial_prob = 1.0 / num_pairs

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
        (phases, haps) = phase(genotypes[i])
        hap_dict = build_haplotype_pair_dictionary(phases)
        gen_dict[genotypes[i]] = hap_dict
        for hap in haps:
            unique_haps.add(hap)

    unique_haps = list(unique_haps)
    num_unique_haps = len(unique_haps)
    initial_hap_prob = 1.0 / num_unique_haps
    for i in range(0, num_unique_haps):
        unique_hap_dict[unique_haps[i]] = initial_hap_prob

    return (gen_dict, unique_hap_dict)

# map all unique haplotypes to their corresponding initial probabilities
def build_unique_haplotype_dictionary(haplotype_dictionary):
    haps = set()
    for (key_gen, value_hap_pair) in haplotype_dictionary:
        hap_pair = value_hap_pair.key()





# def main():
#     genotypes = get_genotypes("/Users/dnicolaou/Desktop/test_data_2.txt")
#     window = genotypes[0][0:51]
#
# main()
