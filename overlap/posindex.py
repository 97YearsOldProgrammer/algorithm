import sys
import json

def read_features(file, size):
	features = {}
	with open(file) as fp:
		for line in fp.readlines():
			chrom, beg, end = line.split()
			beg = int(beg)
			end = int(end)
			bx = beg // size
			ex = end // size

			if chrom not in features: features[chrom] = []
			if len(features[chrom]) <= ex:
				for i in range(len(features[chrom]), ex +1):
					features[chrom].append([])

			for i in range(bx, ex+1):
				features[chrom][i].append( (int(beg), int(end)) )

	return features

def overlaps(f1, f2):
	b1, e1 = f1
	b2, e2 = f2
	if b2 <= e1 and e2 >= b1: return True
	return False

size = 10000 # best value is what?
f1d = read_features(sys.argv[1], size)
f2d = read_features(sys.argv[2], size)

seen = {}
for chrom in f1d:
	for i in range(len(f1d[chrom])):
		for f1 in f1d[chrom][i]:
			if i < len(f2d[chrom]):
				for f2 in f2d[chrom][i]:
					if overlaps(f1, f2):
						b1, e1 = f1
						b2, e2 = f2
						if (chrom, b1, e1, b2, e2) in seen: continue
						print(f'{chrom}:{b1}-{e1} {chrom}:{b2}-{e2}')
						seen[ (chrom, b1, e1, b2, e2) ] = True