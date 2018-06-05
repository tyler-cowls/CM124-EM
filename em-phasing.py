import sys
import time

steady = 0
transient = 1

# given an input genotype string, generate a list of all possible haplotypes
def get_haplotypes(gen):
	haps = []
	# recursively generate haplotypes which are compatible with the genotype
	# initially pass genotype in as tail, and SNP to head after processing it
	# branches into 2 recursive calls at every 1
	def parse(head, tail):
		# base case if tail is empty, finished generating haplotype, so append to list of haplotypes
		if len(tail) == 0:
			haps.append(head);
		else:
			# if genotype's SNP is 0, haplotype's SNP will also be 0
			if tail[0] == '0':
				parse(head+'0', tail[1:])
			elif tail[0] == '2': # if genotype's SNP is 2, haplotype's SNP will be 1
				parse(head+'1', tail[1:])
			else: # if genotype's SNP is 1, branch into two possible haplotypes (one with 1 at SNP, another with 0 at SNP)
				parse(head+'1', tail[1:])
				parse(head+'0', tail[1:])
	parse('', gen)
	return haps

# given a genotype, return all possible phases (pairs of compatible haplotypes)
def phase(gen):
	# get the list of haplotypes associated with genotype gen
	haps = get_haplotypes(gen)
	phases = []
	num_haps = len(haps)
	if num_haps == 1:
		phases.append((haps[0], haps[0]))
	else: # pair the first and last haplotypes together and work inwards
		for i in range(0, (int)(num_haps/2)):
			phases.append((haps[i], haps[-(i+1)]))
	return (phases, haps)

# given a list of genotypes, return a list of the phases, their associated probabilities, 
# and a dictionary of unique haplotypes -> haplotype probabilities
def get_phase_prob_pairs(genotype_list):
	num_gens = len(genotype_list)
	phase_list = []
	phase_probs = []
	unique_haps = set()
	unique_hap_dict = {}
	# for each genotype, get the phases an haplotypes
	for i in range(0, num_gens):
		phases, haps = phase(genotype_list[i])
		phase_list.append(phases)
		probs = []
		num_phases = len(phases)
		for j in range(0, num_phases):
			probs.append(0.0) # initialize phase probability to 0, so can calculate later using haplotype probabilities
		phase_probs.append(probs)
		for hap in haps:
			unique_haps.add(hap)

	unique_haps = list(unique_haps)
	num_unique_haps = len(unique_haps)
	initial_hap_prob = 1.0 / num_unique_haps # assume uniform distribution of haplotype probabilities initially
	# add each unique haplotype to the dictionary with value (steady, transient) probability tuple
	for i in range(0, num_unique_haps):
		unique_hap_dict[unique_haps[i]] = (initial_hap_prob, 0) # (steady, transient) probabilities
		# steady probability persists throughout calculation of next iteration's probability, which is updated in transient

	return (phase_list, phase_probs, unique_hap_dict)

# given a list of genotypes, return a list lists of individuals' phases and their associated probabilities
def EM(genotype_list):
	phase_list, phase_probs, haplotypes = get_phase_prob_pairs(genotype_list)
	num_gens = len(genotype_list)
	iterations = 10 # TODO: update appropriately

	for i in range(0, iterations): # for each iteration
		for g in range(0, num_gens): # for each genotype
			sumTempProb = 0
			num_phases = len(phase_list[g])
			for p in range(0, num_phases): # for each phase, calculate conjunctive probability of its 2 haplotypes
				prob = haplotypes[phase_list[g][p][0]][steady]*haplotypes[phase_list[g][p][1]][steady]
				sumTempProb += prob
				phase_probs[g][p] = prob # temporarily update phase_prob to conjunctive probability (we don't know sumTempProb yet)

			for p in range(0, num_phases):
				prob = phase_probs[g][p] / sumTempProb
				phase_probs[g][p] = prob # re-update phase_prob to actual probability after dividing by sumTempProb
				temp = haplotypes[phase_list[g][p][0]]
				haplotypes[phase_list[g][p][0]] = (temp[steady], temp[transient]+prob) # persist steady probability, update transient probability
				temp = haplotypes[phase_list[g][p][1]]
				haplotypes[phase_list[g][p][1]] = (temp[steady], temp[transient]+prob) # for both haplotypes in the phase
				
		for hap, probs in haplotypes.items():
			haplotypes[hap] = (probs[transient]/(2*num_gens), 0) # update steady probability for next iteration, reset transient to 0

	return (phase_list, phase_probs)


def windows_EM(in_file, out_file, window_size):
	f = open(in_file, 'r')
	input = f.readlines()

	# window_size = 20
	num_snps = len(input)
	num_windows = (int)(num_snps/window_size)

	for i in range(0, num_windows): # for each window (excluding the last with variable length < window_size)
		offset = i*window_size
		sanitized_input = []
		for r in range(0, window_size): # for each row within the window
			sanitized_input.append(input[offset + r].rstrip().split())
		genotype_list = []
		for c in range(0, len(sanitized_input[0])):
			genotype = ''
			for r in range(0, window_size):
				genotype += sanitized_input[r][c]
			genotype_list.append(genotype)

		phase_list, phase_probs = EM(genotype_list)
		
		max_haps = [] # contains the haplotypes with highest phase probability within each genotype across all genotypes within this window
		num_gens = len(phase_list)
		for g in range(0, num_gens): # for each genotype, find the phase with the highest probability
			num_phases = len(phase_list[g])
			maxProb = -1.0
			maxPhase = ('', '')
			for p in range(0, num_phases): # for each phase, compare probability to max
				if phase_probs[g][p] > maxProb:
					maxPhase = phase_list[g][p]
					maxProb = phase_probs[g][p]
			if maxProb >= 0:
				max_haps.extend(maxPhase) # add phase with highest probability to list of phases
		for c in range(0, window_size):
			row = ''
			for r in range(0, 2*num_gens):
				row += max_haps[r][c] + ' '
			if i==0 and c==0:
				with open(out_file, 'w') as f:
					print(row, file=f)
			else:
				with open(out_file, 'a+') as f:
					print(row, file=f)
		print('completed window ' + str(i) + ' out of ' + str(num_windows))

	# run EM again on final window of variable length which was excluded earlier
	if num_snps%window_size > 0:
		sanitized_input = []
		for r in range(num_windows*window_size, num_snps):
			sanitized_input.append(input[r].rstrip().split())
		genotype_list = []
		for c in range(0, len(sanitized_input[0])):
			genotype = ''
			for r in range(0, num_snps%window_size):
				genotype += sanitized_input[r][c]
			genotype_list.append(genotype)

		phase_list, phase_probs = EM(genotype_list)

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
		for c in range(0, num_snps%window_size):
			row = ''
			for r in range(0, 2*num_gens):
				row += max_haps[r][c] + ' '
			with open(out_file, 'a+') as f:
				print(row, file=f)
	print('completed all windows')


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