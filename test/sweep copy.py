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
        events.append({
            'position': beg,
            'priority': 0,
            'scope': (beg, end),
            'name': (idx, beg, end),
            'overlap': []
        })
        events.append({
            'position': end,
            'priority': 1,
            'name': (idx, beg, end)
        })

if arg.gff2:
    with open(arg.gff2, 'r') as fp:
        for line in fp:
            idx, beg, end = line.strip().split('\t')
            beg = int(beg)
            end = int(end)
            events.append({
                'position': beg,
                'priority': 0,
                'name': (idx, beg, end),
            })
            events.append({
                'position': end,
                'priority': 1,
                'name': (idx, beg, end)
            })

events.sort(key=lambda x: (x['position'], x['priority']))
# Algorithm
overlaps = []
for event in events:
    priority = event['priority']
    name     = event['name']
    if priority == 0:
        overlaps.append(name)
    elif priority == 1:
        overlaps.remove(name)
        name1, beg1, end1 = name
        if overlaps:
            print(f'{name1}\t{beg1}\t{end1}')
            for name2, beg2, end2 in overlaps:
                print(f'\t{name2}\t{beg2}\t{end2}')
            print()
        elif not overlaps: 
            print(f'There is no overlap for {name1}\t{beg1}\t{end1}')
            print()