import gzip
import argparse
import math

def freq(counts, k):
	assert(len(counts) == 4**k) 
	total = 0
	for kmer in counts: total += counts[kmer]
	freq = {}
	for kmer in counts: freq[kmer] = counts[kmer]/total
	return freq
	

parser = argparse.ArgumentParser(
	description='generic ime trainer')
parser.add_argument('master', type=str, metavar='<file>',
	help='path to *.master.txt.gz')
parser.add_argument('--k', type=int, metavar='<int>', default=5,
	required=False, help='k-mer size [%(default)i]')
parser.add_argument('--d', type=int, metavar='<int>', default=5,
	required=False, help='donor length [%(default)i]')
parser.add_argument('--a', type=int, metavar='<int>', default=10,
	required=False, help='acceptor length [%(default)i]')
parser.add_argument('--prox', type=int, metavar='<int>', default=400,
	required=False, help='proximal cutoff [%(default)i]')
parser.add_argument('--dist', type=int, metavar='<int>', default=400,
	required=False, help='distal cutoff [%(default)i]')
parser.add_argument('--maxlen', type=int, metavar='<int>', default=1000,
	required=False, help='maximum intron length [%(default)i]')
parser.add_argument('--minexp', type=float, metavar='<float>', default=100,
	required=False, help='minimum expression [%(default)i]')
arg = parser.parse_args()

p = {} #ime k-mers
q = {} #non ime k-mers
with gzip.open(arg.master, "rt") as fp:
	for line in fp.readlines(): 
		name, beg, end, strand, exp, seq = line.rstrip().split()
		beg = int(beg)
		end = int(end)
		exp = float(exp)
		if end - beg >= arg.maxlen: continue
		if beg < arg.prox and exp > arg.minexp: c = p
		else: c = q
		for i in range(arg.d, len(seq) - arg.a - arg.k + 1):
			kmer = seq[i:i+arg.k]
			if kmer not in c: c[kmer] = 0
			c[kmer] += 1
p = freq(p, arg.k)
q = freq(q, arg.k)
for kmer in p: print(kmer, math.log2(p[kmer]/q[kmer]))
