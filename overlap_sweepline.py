import argparse
import gzip
import optool
from datetime import datetime

parser = argparse.ArgumentParser(description='Overlap finder for transcripts exons/introns')
parser.add_argument('gff', type=str, help='GFF file (can be gzipped).')
arg = parser.parse_args()

start = datetime.now()

feature = ['exon', 'intron']
genes = {}

with gzip.open(arg.gff, 'rt') as fp:
    for line in fp:
        line = line.strip()
        x    = line.split('\t')
        src = x[1]
        if src != 'WormBase': continue
        typ = x[2]
        if typ not in feature: continue
        chm = x[0]
        beg = x[3]
        end = x[4]
        att = x[8]
        if ';' in att:
            att = att.split(';')
            att = att[0]
        att = att.split(':')
        name = att[1]
        if name not in genes: genes[name] = []
        genes[name].append({
            'type': typ,
            'beg': int(beg),
            'end': int(end)
            })

# create events for algorithm
events = [] 
for trxn in genes:
    for feature in genes[trxn]:
        events.append({
            'position': feature['beg'],
            'priority': 0, 
            'scope': (feature['beg'], feature['end']),
            'type': feature['type'],
            'name': trxn
        })
        events.append({
            'position': feature['end'],
            'priority': 1, 
            'scope': (feature['beg'], feature['end']),
            'type': feature['type'],
            'name': trxn
        })
events.sort(key=lambda x: ( x['position'], x['priority'] ) )

# this gonna be a complex data strcuture
overlaps = {} 
# key   ( trxn_name, typ, beg, end )
# value ( trxn_name, typ, beg, end )
# value is a set(make sure no depulicate) with what are overlapped
overlap  = set()
# there is how we know there are overlap
alternative_splicing = {}
# key   (trxn name)
# value (type we have)

# algorithm start
for event in events:
    # unpackle those values
    position = event['position']
    priority = event['priority']
    name     = event['name']
    scope    = event['scope']
    typ      = event['type']
    beg, end = scope
    # this is how we are gonna know component of overall structure
    if name not in alternative_splicing: alternative_splicing[name] = set()
    # this is what would be value 
    index    = ( name, typ, beg, end )
    # so that we can know how different introns and exons overlap with each other
    if index not in overlaps: overlaps[index] = set()
    # start the algorithm
    if priority == 0:
        # append overlap 
        for olp in overlap:
            overlaps[olp].add( ( name, typ, beg, end ))
            overlaps[index].add( olp )  # Also add existing features to overlaps of current feature
        # add overlap index
        overlap.add(index)
        # how we know component of different transcription
        alternative_splicing[name].add( ( typ, beg, end ) )
    elif priority == 1:
        overlap.remove(index)

output = {}
# n as name, t as type, b as beg, e as end
for tuple1 in overlaps:
    n1, t1, b1, e1 = tuple1
    if n1 not in output: output[n1] = []
    for tuple2 in overlaps[tuple1]:
        n2, t2, b2, e2 = tuple2
        output[n1].append( (t1, b1, e1, t2, n2, b2, e2) )

# list for no those don't have output
no_overlap = []

for transcripts in output:
    # sorting output if there do have overlap
    if output[transcripts]:
        output[transcripts].sort(key=lambda x: (x[1],x[2],x[5],x[6]) )
        print(f'{transcripts}')
        # final output
        for overlap in output[transcripts]:
            t1, b1, e1, t2, n2, b2, e2 = overlap
            print(f'\t{t1:<10}\t{b1:<10}\t{e1:<10}\t{t2:<10}\t{n2:<15}\t{b2:<10}\t{e2:<10}')
    # label those who don't have overlap
    else: no_overlap.append(transcripts)

for transcript in no_overlap:
    print(f'{transcript}\tThere is no overlap')

end = datetime.now()
print(f"Total processing time: {end - start}")
