import gzip
import sys

#thresholds: top 100 ime (>464.7)
# min and max motif widths (-minw, -maxw), max number of motifs (-nmotifs)

def read_model(filename):
	model = {}

	fp = None
	if   filename.endswith('.gz'): fp = gzip.open(filename, 'rt')
	else:                          fp = open(filename)

	for line in fp.readlines():
		if line.startswith('#'): continue
		kmer, imeter = line.split()
		model[kmer] = float(imeter)

	return model

def score_intron(model, seq, k, d=5, a=10):
	s = 0
	for i in range(d, len(seq) -k + 1 -a):
		kmer = seq[i:i+k]
		if kmer in model: s += model[kmer]
	return s
	
ime_model = read_model("../datacore/project_imeter/dm_ime_model")

min_ime = 464.7
keep    = {}

with gzip.open(sys.argv[1], "rt") as fp:
	for line in fp.readlines():
		line = line.rstrip()
		id, b, e, s, x, seq = line.split()
		b = int(b)
		e = int(e)
		ime = score_intron(ime_model, seq, 4)
		if ime < min_ime: continue
		if seq not in keep: keep[seq] = f">{id}.{b}.{e}"

for seq in keep: 
	print(keep[seq])
	print(seq)
