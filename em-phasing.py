import sys
import time

steady = 0
transient = 1

def get_haplotypes(gen):
	haps = []
	def parse(head, tail):
		if len(tail) == 0:
			haps.append(head);
		else:
			if tail[0] == '0':
				parse(head+'0', tail[1:])
			elif tail[0] == '2':
				parse(head+'1', tail[1:])
			else:
				parse(head+'1', tail[1:])
				parse(head+'0', tail[1:])
	parse('', gen)
	return haps

def phase(gen):
	haps = get_haplotypes(gen)
	phases = []
	num_haps = len(haps)
	if num_haps == 1:
		phases.append((haps[0], haps[0]))
	else:
		for i in range(0, (int)(num_haps/2)):
			phases.append((haps[i], haps[-(i+1)]))
	return (phases, haps)

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
		phases, haps = phase(genotype_list[i])
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
	iterations = 15 # TODO: update appropriately

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


def windows_EM(in_file, out_file, window_size):
	f = open(in_file, 'r')
	input = f.readlines()

	# window_size = 20
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

		phase_list, phase_probs, haplotypes = EM(genotype_list)
		
		max_haps = []
		num_gens = len(phase_list)
		for g in range(0, num_gens): # for each genotype
			num_phases = len(phase_list[g])
			maxProb = -1.0
			maxPhase = ('', '')
			for p in range(0, num_phases): # for each phase
				if phase_probs[g][p] > maxProb:
					maxPhase = phase_list[g][p]
					maxProb = phase_probs[g][p]
			if maxProb >= 0:
				max_haps.extend(maxPhase)
		total_max_haps.append(max_haps)
		print('completed window ' + str(i) + ' out of ' + str(num_windows))

	for r in range(num_windows*window_size, num_snps):
		input[r] = input[r].rstrip().split()
	genotype_list = []
	for c in range(0, len(input[0])):
		genotype = ''
		for r in range(num_windows*window_size, num_snps):
			genotype += input[r][c]
		genotype_list.append(genotype)

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

	final_max_haps = total_max_haps[0]
	l = len(total_max_haps)
	l2 = len(final_max_haps)
	for i in range(1, l):
		for j in range(0, l2):
			final_max_haps[j] += total_max_haps[i][j]

	l3 = len(final_max_haps[0])
	for i in range(0, l3):
		row = ''
		for h in final_max_haps:
			row += h[i] + ' '
		if i==0:
			with open(out_file, 'w') as f:
				print(row, file=f)
		else:
			with open(out_file, 'a+') as f:
				print(row, file=f)


def main():
	in_file = ''
	out_file = ''
	window_size = 15 # default
	if len(sys.argv) == 4:
		in_file = sys.argv[1]
		out_file = sys.argv[2]
		window_size = (int)(sys.argv[3])
	elif len(sys.argv) == 3:
		in_file = sys.argv[1]
		out_file = sys.argv[2]
	elif len(sys.argv) == 2:
		in_file = sys.argv[1]
		out_file = 'out.txt'
	else:
		print('Error: incorrect number of arguments ')
		return -1

	start = time.time()
	windows_EM(in_file, out_file, window_size)
	print('Total time: ', time.time() - start)
	

if __name__ == "__main__":
    main()