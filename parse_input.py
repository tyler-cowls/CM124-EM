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
def build_haplotype_dictionary(haplotype_pair_list):
    num_pairs = len(haplotype_pair_list)
    haplotype_dict = {}
    initial_prob = 1.0 / num_pairs
    print initial_prob

    for i in range(0, num_pairs):
        haplotype_dict[haplotype_pair_list[i]] = initial_prob

    return haplotype_dict

def main():
    genotypes = get_genotypes("/Users/dnicolaou/Desktop/test_data_2.txt")
    window = genotypes[0][0:51]

main()
