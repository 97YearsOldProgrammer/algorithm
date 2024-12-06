import argparse

parser = argparse.ArgumentParser(description='Overlap gene reporter')
parser.add_argument('gff1', type=str, help='GFF1 file')
parser.add_argument('--gff2', type=str, help='GFF2 file')
arg = parser.parse_args()

events = []
with open(arg.gff1, 'r') as fp:
    for line in fp:
        idx, beg, end = line.strip().split('\t')
        beg = int(beg)
        end = int(end)
        events.append( (idx, beg, end, 0) )
        events.append( (idx, beg, end, 1) )

if arg.gff2:
    with open(arg.gff2, 'r') as fp:
        for line in fp:
            idx, beg, end = line.strip().split('\t')
            beg = int(beg)
            end = int(end)
            events.append( (idx, beg, end, 0) )
            events.append( (idx, beg, end, 1) )

events.sort(key=lambda x: ( x[1], x[3] ))
overlaps = []
for idx, beg, end, typ in events:
    output = (idx, beg, end)
    if typ == 0:
        overlaps.append(output)
    elif typ == 1:
        overlaps.remove(output)
        if overlaps:
            print(f'{idx}\t{beg}\t{end}')
            for name, beg1, end1 in overlaps:
                print(f'\t{name}\t{beg1}\t{end1}')
            print()
        elif not overlaps: 
            print(f'There is no overlap for {idx}\t{beg}\t{end}')
            print()

