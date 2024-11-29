import gzip
import sys

# what we would use in gff_reader for finding name of gene
# this is important for know combination of exons and introns
# ex: attribute : eight column of gff
# key what we want to find stuff behind ID; key=ID

def gene_namer(attributes, key):
    attributes = attributes.strip().split(';')
    for attribute in attributes:
        attribute = attribute.strip()
        if attribute.startswith(f'{key}='):
            name = attribute.split('=', 1)[1]
            return name
    return ''

# easy version of overlap to see overlap gene 
# input: gff file
def overlap_gene_gff(filename):
    """ credit to Dr.Korf """
    if       filename == '-':                  fp = sys.stdin
    elif     filename.endswith('.gz'):         fp = gzip.open(filename, 'rt')
    else:                                      fp = open(filename)
    gff  = {}
    while True:
        line = fp.readline()
        if line == '':break
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
    fp.close()

'''
# extract information for sweep line algorithm from gff
def gff_reader(filename):
    """ credit to Dr.Korf """
    if       filename == '-':                  fp = sys.stdin
    elif     filename.endswith('.gz'):         fp = gzip.open(filename, 'rt')
    else:                                      fp = open(filename)
    gff  = []
    set  = set() 
    gene_type = ['mRNA', 'piRNA','snoRNA', 
                 'miRNA','pre_miRNA',
                 'pseudogenic_transcript'
                 'tRNA','pseudogenic_tRNA',
                 'ncRNA', 'circular_ncRNA' , 'nc_primary_transcript', ]
    forward_collector = False     # so that we can know whenever we need to collect proper information
    reverse_collector = False     # we split them into forward and reverse strand
    
    while True:
        line = fp.readline()
        if line == '':break
        x    = gff_line_splitter(line)
        if x['chm'] not in gff: gff[ x['chm'] ] = []
        """ we only need line from WormBase"""
        if x['src'] == 'WormBase':
            """ the actual line below gene of WormBase contain everything we need"""
            if x['typ'] == 'gene':
                if   x['pol'] == '+':
                    forward_strand = []
                    forward_collector = True
                elif x['pol'] == '-':
                    reverse_strand = []
                    reverse_collector = True
            """ where we start collecting specific ID for different transcript under gene"""
            if   forward_collector:
                if x['typ'] in gene_type: 
                    name = gene_namer(x['inf'])
                    forward_strand.append(name)
                
            elif reverse_collector:
                if x['typ'] in gene_type: 
                    name = gene_namer(x['inf'])
                    reverse_strand.append(name)

            if forward_strand:
                if id in forward_strand in line:
                    gff[chm].append( ())
            
            elif reverse_strand:


                    beg = x[3]
                    end = x[4]
                    pol = x[6]
                    inf = x[8]
                # get ID of gene
                id  = gene_namer(inf)
            gff[chm].append( (int(beg), int(end), typ) )
        else: continue
'''