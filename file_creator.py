import argparse
import random

parser = argparse.ArgumentParser(description='Overlap Test File Creator')
parser.add_argument('output_file', type=str, help='Output filename')
parser.add_argument('--olp', type=int, default=5,
    help='number of overlaps [%(default)i]')
parser.add_argument('--noolp', type=int, default=1,
    help='number of nonoverlaps [%(default)i]')
arg = parser.parse_args()

with open(arg.output_file, 'w') as outfile:
    for i in range(arg.olp):
        name = f'Trxn{i+1}'
        beg_idx = None
        for j in range(3):
            if not beg_idx:
                beg  = random.randint(1, 1000)
                end  = random.randint(beg, 1000)
                beg_idx = end
            elif beg_idx:
                beg  = random.randint(beg_idx, 1000)
                end  = random.randint(beg, 1000)
                beg_idx = end
            if beg + 20 > 1000: break
            outfile.write(f'{name}\t{beg}\t{end}\n')

    for i in range(arg.noolp):
        name = f'Trxn{arg.olp + i + 1}'
        beg = random.randint(10000, 50000)
        end = random.randint(beg, 50000)
        outfile.write(f'{name}\t{beg}\t{end}\n')

