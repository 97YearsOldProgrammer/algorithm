import argparse
import optool

parser = argparse.ArgumentParser(description='overlap gene reporter')
parser.add_argument('gff', type=str, help='GFF file (gzip file)')
arg = parser.parse_args()

fp = open(arg.gff)

gene = []
gene_type = ['mRNA', 'piRNA','snoRNA', 
                 'miRNA','pre_miRNA',
                 'pseudogenic_transcript'
                 'tRNA','pseudogenic_tRNA',
                 'ncRNA', 'circular_ncRNA' , 'nc_primary_transcript', ]
transcripts = {}

gene_reader = False
while True:
    line = fp.readline()
    if line == '':break
    x = line.split()
    typ = x[2]
    att = x[8]

    if gene_reader:

        if typ in gene_type:
            transcript = optool.gene_namer(att, 'Name')
            if transcript not in transcripts[current_gene]:
                transcripts[current_gene][transcript] = []

        if typ == 'gene':
            gene_reader = False

    if 'Y74C9A.2' in att:
        current_gene = 'Y74C9A.2'
        if current_gene not in transcripts:
            transcripts[current_gene] = {}
        gene_reader = True
    

fp.close()
print(transcripts)

'''
    x = line.split()
    chm = x[0]
    if chm not in gff: gff[chm] = []
    typ = x[2]
    beg = x[3]
    end = x[4]
    pol = x[6]
    att = x[8]
    if typ == 'gene': 
        name = gene_namer(att, 'sequence_name')
        gff[chm].append( (beg, end, pol, name) )
    return gff
'''
