from parse_input import *

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



def main():
	gens = ['21110', '10001', '22221']
	#gens = ['202020']

	# print(phase(gens[2]))
	# print(build_haplotype_pair_dictionary(phase(gens[2])))
	(dict, haps) = build_dictionaries(gens)
	print(dict)
	print(haps)

if __name__ == "__main__":
	main()
