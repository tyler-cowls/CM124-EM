import create_phases
import sys

steady = 0
transient = 1

# outputs list of genotypes
def get_genotypes(file_path):
	f = open(file_path, 'r')
	snps = f.readlines()
    
	for i in range(0, len(snps)):
		snps[i] = snps[i].rstrip()
		snps[i] = snps[i].split()

	genotype_list = list()

	for i in range(0, len(snps[0])):
		genotype = ""
        
		for j in range(0, len(snps)):
			genotype += snps[j][i]

		genotype_list.append(genotype)

	return genotype_list

def get_phase_prob_pairs(genotype_list):
	num_gens = len(genotype_list)
	phase_list = []
	phase_probs = []
	unique_haps = set()
	unique_hap_dict = {}
	for i in range(0, num_gens):
		phases, haps = create_phases.phase(genotype_list[i])
		phase_list.append(phases)
		probs = []
		num_phases = len(phases)
		for j in range(0, num_phases):
			probs.append(0.0) #1.0/num_phases
		phase_probs.append(probs)
		for hap in haps:
			unique_haps.add(hap)

	unique_haps = list(unique_haps)
	num_unique_haps = len(unique_haps)
	initial_hap_prob = 1.0 / num_unique_haps
	for i in range(0, num_unique_haps):
		unique_hap_dict[unique_haps[i]] = (initial_hap_prob, 0) # (steady, transient) probabilities

	return (phase_list, phase_probs, unique_hap_dict)

def EM(genotype_list):
	phase_list, phase_probs, haplotypes = get_phase_prob_pairs(genotype_list)
	num_gens = len(genotype_list)
	iterations = 50 # TODO: update appropriately

	for i in range(0, iterations): # for each iteration
		for g in range(0, num_gens): # for each genotype
			sumTempProb = 0
			num_phases = len(phase_list[g])
			for p in range(0, num_phases): # for each phase, calculate conjunctive probability of its 2 haplotypes
				prob = haplotypes[phase_list[g][p][0]][steady]*haplotypes[phase_list[g][p][1]][steady]
				sumTempProb += prob
				phase_probs[g][p] = prob

			for p in range(0, num_phases):
				prob = phase_probs[g][p] / sumTempProb
				phase_probs[g][p] = prob
				temp = haplotypes[phase_list[g][p][0]]
				haplotypes[phase_list[g][p][0]] = (temp[steady], temp[transient]+prob)
				temp = haplotypes[phase_list[g][p][1]]
				haplotypes[phase_list[g][p][1]] = (temp[steady], temp[transient]+prob)
				
		for hap, probs in haplotypes.items():
			haplotypes[hap] = (probs[transient]/(2*num_gens), 0)

	return (phase_list, phase_probs, haplotypes)

def non_windows_EM(file_path):
	genotype_list = get_genotypes(sys.argv[1])

	phase_list, phase_probs, haplotypes = EM(genotype_list)

	final_haps = []
	num_gens = len(phase_list)
	for i in range(0, num_gens): # for each genotype
		num_phases = len(phase_list[i])
		maxProb = -1.0
		maxPhase = ('', '')
		for j in range(0, num_phases): # for each phase
			if phase_probs[i][j] > maxProb:
				maxPhase = phase_list[i][j]
				maxProb = phase_probs[i][j]
		if maxProb >= 0:
			final_haps.extend(maxPhase)

	gen_len = len(genotype_list[0])
	for i in range(0, gen_len):
		row = ''
		for h in final_haps:
			row += h[i] + ' '
		print(row)

def partial_windows_EM(file_path):
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
		phase_list, phase_probs, haplotypes = EM(genotype_list)
		
		max_haps = []
		num_gens = len(phase_list)
		for i in range(0, num_gens): # for each genotype
			num_phases = len(phase_list[i])
			maxProb = -1.0
			maxPhase = ('', '')
			for j in range(0, num_phases): # for each phase
				if phase_probs[i][j] > maxProb:
					maxPhase = phase_list[i][j]
					maxProb = phase_probs[i][j]
			if maxProb >= 0:
				max_haps.extend(maxPhase)
		total_max_haps.append(max_haps)

	print(total_max_haps)

def windows_EM(file_path):
	f = open(file_path, 'r')
	input = f.readlines()

	window_size = 50
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

		#print(genotype_list)
		phase_list, phase_probs, haplotypes = EM(genotype_list)
		
		max_haps = []
		num_gens = len(phase_list)
		for i in range(0, num_gens): # for each genotype
			num_phases = len(phase_list[i])
			maxProb = -1.0
			maxPhase = ('', '')
			for j in range(0, num_phases): # for each phase
				if phase_probs[i][j] > maxProb:
					maxPhase = phase_list[i][j]
					maxProb = phase_probs[i][j]
			if maxProb >= 0:
				max_haps.extend(maxPhase)
		total_max_haps.append(max_haps)

	for r in range(num_windows*window_size, num_snps):
		input[r] = input[r].rstrip().split()
	genotype_list = []
	for c in range(0, len(input[0])):
		genotype = ''
		for r in range(num_windows*window_size, num_snps):
			genotype += input[r][c]
		genotype_list.append(genotype)

	#print(genotype_list)
	phase_list, phase_probs, haplotypes = EM(genotype_list)

	max_haps = []
	num_gens = len(phase_list)
	for i in range(0, num_gens): # for each genotype
		num_phases = len(phase_list[i])
		maxProb = -1.0
		maxPhase = ('', '')
		for j in range(0, num_phases): # for each phase
			if phase_probs[i][j] > maxProb:
				maxPhase = phase_list[i][j]
				maxProb = phase_probs[i][j]
		if maxProb >= 0:
			max_haps.extend(maxPhase)
	total_max_haps.append(max_haps)

	#print(total_max_haps)

	final_max_haps = total_max_haps[0]
	l = len(total_max_haps)
	l2 = len(final_max_haps)
	for i in range(1, l):
		for j in range(0, l2):
			final_max_haps[j] += total_max_haps[i][j]

	#print(final_max_haps)

	with open('out.txt', 'w') as f:
		print('', file=f)
	l3 = len(final_max_haps[0])
	for i in range(0, l3):
		row = ''
		for h in final_max_haps:
			row += h[i] + ' '
		if i==0:
			with open('out.txt', 'w') as f:
				print(row, file=f)
		else:
			with open('out.txt', 'a+') as f:
				print(row, file=f)


def main():
	#non_windows_EM(sys.argv[1])
	windows_EM(sys.argv[1])
	

if __name__ == "__main__":
    main()