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

def main():
    genotypes = get_genotypes("/Users/dnicolaou/Desktop/test_data_2.txt")
    window = genotypes[0][0:51]
    print generate_haplotypes(window)

main()
