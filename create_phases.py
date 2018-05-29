from parse_input import build_haplotype_dictionary

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

def phase(gen):
	parse('', gen)
	phases = []
	for i in range(0, (int)(len(haps)/2)):
		phases.append((haps[i], haps[-(i+1)]))
	return phases



def main():
	gens = ['21110', '10001', '22221']

	# print(phase(gens[0]))
	print(build_haplotype_dictionary(phase(gens[0])))

if __name__ == "__main__":
	main()
