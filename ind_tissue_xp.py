import gzip
import sys
import math
import statistics

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

#get raw counts
xp_raw = []

for i in range(0,11):
	xp_raw.append([])

with gzip.open(sys.argv[1], "rt") as fp:
	for line in fp.readlines():
		line = line.rstrip()
		f = line.split()
		for i in range(11):
			xp_raw[i].append(int(f[i+4]))

#get the avg
xp_avg = []	
for i in range(len(xp_raw)):
	xp_avg.append(statistics.mean(xp_raw[i]))

print(xp_avg)

all = []
tsp = []
rst = []

#create lists of tissue specificity
with gzip.open(sys.argv[1], "rt") as fp:
	for line in fp.readlines():
		line = line.rstrip()
		f = line.split()
		if int(f[1]) > 400: continue
		tissues = []
		for i in range(11):
			tissues.append(int(f[i+4])/xp_avg[i])
		h = shape(tissues)
		if h == None: continue
		if h < 1: tsp.append(f[-1])
		elif h < 2: rst.append(f[-1])
		else: all.append(f[-1])
		
		
print(len(all), len(rst), len(tsp))

		
		


	


