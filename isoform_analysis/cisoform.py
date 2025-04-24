import os
import glob
import argparse

import isoform2
import hints

parser = argparse.ArgumentParser(
	description='Alternative isoform generator')
parser.add_argument('dir', type=str, metavar='<fastas dir>',
	help='dir for the small gene set')
parser.add_argument('--min_intron', required=False, type=int, default=35,
	metavar='<int>', help='minimum length of intron [%(default)i]')
parser.add_argument('--min_exon', required=False, type=int, default=25,
	metavar='<int>', help='minimum length exon [%(default)i]')
parser.add_argument('--flank', required=False, type=int, default=99,
	metavar='<int>', help='genomic flank on each side [%(default)i]')
parser.add_argument('--limit_fasta', required=False, type=int, default=1,
	metavar='<int>', help='limit number of fasta [%(default)i]')
arg = parser.parse_args()

def run(dir):
    fastas = glob.glob(os.path.join(dir, "*.fa"))
    fcount  = 0
    
    for fasta in fastas:
        fcount += 1
        if fcount > arg.limit_fasta: break
        
        with open(fasta, 'r') as file:
            name, seq  = next(isoform2.read_fasta(fasta))
            dons, accs = isoform2.gtag_sites(seq, arg.flank, arg.min_exon)
            count      = hints.countiso(dons, accs, arg.min_intron, arg.min_exon)
            isoform, info = hints.all_possible2(seq, arg.min_intron, arg.min_exon, 30, arg.flank)
            print(f'{name}\t{count}\t{len(isoform)}')
            
run(arg.dir)
