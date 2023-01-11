import math
import sys
import statistics
import gzip

'''
compare low, medium, and high expression levels in tissue specific, restricted, and all specificity levels
'''

def shape(vals):
	s = sum(vals)
	if s == 0:
		return None
	p = []
	for val in vals:
		p.append(val/s)
	return entropy(p)

def entropy(ps):
	h = 0
	for p in ps:
		if p != 0:
			h -= p*math.log2(p)
	return h

def read_model(filename):
	model = {}

	fp = None
	if   filename.endswith('.gz'): fp = gzip.open(filename, 'rt')
	else:                          fp = open(filename)

	for line in fp.readlines():
		if line.startswith('#'): continue
		#kmer, imeter, prox, dist = line.split()
		kmer, imeter = line.split()
		model[kmer] = float(imeter)

	return model

def score_intron(model, seq, k, d=5, a=10):
	s = 0
	for i in range(d, len(seq) -k + 1 -a):
		kmer = seq[i:i+k]
		if kmer in model: s += model[kmer]
	return s
	
	
def categorize(bounds, val):
	for i,b in enumerate(bounds):
		if val < b: return i
	return len(bounds) 
	
ime_model = read_model("../../imeter/v1/imeter1.model")
#ime_model = read_model("dm_ime_model")

#get raw counts
xp_raw = []

for i in range(11): xp_raw.append([])

ibs  = [200, 400, 800, 1600]
ds   = [100, 200, 400, 800]
imes = [-10, 0, 10, 20, 40]

with gzip.open(sys.argv[1], "rt") as fp:
	for line in fp.readlines():
		line = line.rstrip()
		f = line.split()
		#if int(f[1]) > 400: continue #proximal introns
		tissues = []
		for i in range(11):
			tissues.append(float(f[i+4]))
		x = sum(tissues)/208622 #max_xp
		if x == 0: continue
		x = math.log2(x)
		h = shape(tissues)
		ib = categorize(ibs,int(f[1]))
		ie = int(f[2])
		d = categorize(ds, (int(f[2]) - int(f[1])))
		ime = categorize(imes, score_intron(ime_model, f[-1], 5))
		#ib = math.log2(ib)
		print(f[0], h, x, ib, d, ime)
		
		
'''
categories:
	ib: (200, 400, 800, 1600)
	d: (100, 200, 400, 800)
	ime: (-10, 0, 10, 20, 40)
'''

		


