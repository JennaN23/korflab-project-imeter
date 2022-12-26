import gzip
import sys

#thresholds: length, distance, and expression 
# min and max motif widths (-minw, -maxw), max number of motifs (-nmotifs)

minlen = 200
maxlen = 1000
maxdist = 500
minx = 50000  
keep = {}

with gzip.open(sys.argv[1], "rt") as fp:
	for line in fp.readlines():
		line = line.rstrip()
		id, b, e, s, x, seq = line.split()
		b = int(b)
		e = int(e)
		x = float(x)
		if b > maxdist: continue
		if x < minx: continue
		if e-b < minlen: continue
		if e-b > maxlen: continue
		if seq not in keep: keep[seq] = f">{id}.{b}.{e}"

for seq in keep: 
	print(keep[seq])
	print(seq)
