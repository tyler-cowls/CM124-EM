# functions to parse input file into dictionaries

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

def build_snp_list(num_haplotypes, haplotype_list):

    if num_haplotypes < 1:
        return haplotype_list

    for i in range(0, num_haplotypes/2):
        haplotype_list[i] += 1
    for i in range((num_haplotypes/2) + 1, num_haplotypes):
        haplotype_list[i] += 0

    build_snp_list(num_haplotypes/2, haplotype_list[0:(num_haplotypes/2])+1)
    build_snp_list(num_haplotypes/2, haplotype_list[num_haplotypes/2:num_haplotypes])


# generate all haplotypes corresponding to a given genotype_window
def generate_haplotypes(genotype_window):
    unknown_positions = list()

    # get list of all unknown positions in genotype
    for i in range(0, len(genotype)):
        if genotype[i] == "1":
            unknown_positions.append(i)

    num_unknown_positions = len(unknown_positions)
    num_haplotypes = 2^num_unknown_positions




    return len(unknown_positions)

def main():
    genotypes = get_genotypes("/Users/dnicolaou/Desktop/test_data_2.txt")
    window = genotypes[0][0:51]
    print generate_haplotypes(window)

main()
