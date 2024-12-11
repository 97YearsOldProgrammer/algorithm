import argparse
import random
import sys

parser = argparse.ArgumentParser(
	description='Create random sequence features')
parser.add_argument('--features', type=int, metavar='<int>', required=False,
	default=10, help = 'number of features to generate [%(default)i]')
parser.add_argument('--chroms', type=int, metavar='<int>', required=False,
	default=20, help = 'number of chromosomes [%(default)i]')
parser.add_argument('--length', type=int, metavar='<int>', required=False,
	default=1000, help = 'lengths of chromosomes [%(default)i]')
parser.add_argument('--size', type=int, metavar='<int>', required=False,
	default=500, help = 'average size of feature [%(default)i]')
arg = parser.parse_args()

for i in range(arg.features):
	chrom = random.randint(1, arg.chroms)
	beg = random.randint(1, arg.length)
	end = beg + random.randint(1, arg.size)
	print(f'chr{chrom}\t{beg}\t{end}')