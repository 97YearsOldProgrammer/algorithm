import sys

def read_features(file):
	features = {}
	with open(file) as fp:
		for line in fp.readlines():
			chrom, beg, end = line.split()
			if chrom not in features: features[chrom] = []
			features[chrom].append( (int(beg), int(end)) )
	return features

def overlaps(f1, f2):
	b1, e1 = f1
	b2, e2 = f2
	if b2 <= e1 and e2 >= b1: return True
	return False

f1d = read_features(sys.argv[1])
f2d = read_features(sys.argv[2])

for chrom in f1d:
	for f1 in f1d[chrom]:
		for f2 in f2d[chrom]:
			if overlaps(f1, f2):
				b1, e1 = f1
				b2, e2 = f2
				print(f'{chrom}:{b1}-{e1} {chrom}:{b2}-{e2}')