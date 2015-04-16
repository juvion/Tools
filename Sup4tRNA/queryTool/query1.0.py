#!/usr/bin/python
import sys,re,itertools

# index = {"-1":1, }

def query(seq, loc, base):
	base = base.upper()
	base_list = list(base)
	if base_list[-1] == "U":
		base_list[-1] = "T"
		base = ''.join(base_list)
	if base[0] != "!":
		if seq[int(loc)-1] == base[-1]:
			return 1
		else:
			return 0
	elif base[0] == "!":
		if seq[int(loc)-1] == base[-1]:
			return 0
		else:
			return 1


if len(sys.argv) == 3:
	#open a output file for all filtering lines
	fo_all = open(sys.argv[2] + "_all"  + ".out", 'wa')
	#read in filtering lines and creat output file for each filtering line
	with open(sys.argv[1]) as fi:
		i = 0
		loc = {}
		base = {}
		fo = {}
		for line in fi:
			line = line.lstrip()
			if line[0] == "#":
				continue
			else:
				words = line.split()
				loc[i] = int(words[0])
				base[i] = words[1]
				fo[i] = open(sys.argv[2] + "_" + str(i+1) + ".out", 'wa')
				i += 1

	with open("tRNAs.mfa") as fi:
		for line1, line2 in itertools.izip_longest(*[fi]*2):
			all_q = []
			seq = line2.strip()
			for i in range(len(loc)):
				single_q = query(seq, loc[i], base[i])
				all_q.append(single_q)
				if single_q:
					fo[i].write(line1 + line2)
			if all(all_q):
				fo_all.write(line1 + line2)

	fo_all.close()
	for i in range(len(fo)):
		fo[i].close()
else:
	print "Error: can\'t find file or read data. Usage: query1.0.py file_input prefix_output"
		
