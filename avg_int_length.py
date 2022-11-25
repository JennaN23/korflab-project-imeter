import gzip
import sys
import math


total = 0
count = 0
lengths = [] 
nt = 0
a = 0
c = 0
g = 0
t = 0
xp = [0]*30
max_xp = 0

with gzip.open(sys.argv[1], "rt") as fp:
	for line in fp.readlines():
		line = line.rstrip()
		name, ib, ie, strand, score, seq = line.split()
		score = float(score)
		if score > max_xp: max_xp = score
		if score != 0:
			xp[int(math.log2(score))] += 1
		count += 1
		total += len(seq)
		lengths.append(len(seq))
		a += seq.count("A")
		c += seq.count("C")
		g += seq.count("G")
		t += seq.count("T")
		nt += len(seq)
lengths.sort()
print(max_xp)
print(lengths[int(len(lengths)/2)])
print(total/count)
print(a/nt)
print(c/nt)
print(g/nt)
print(t/nt)
for i in range(len(xp)):
	print(2**i, xp[i])

		
