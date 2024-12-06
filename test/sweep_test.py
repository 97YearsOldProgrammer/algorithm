import argparse
from datetime import datetime

parser = argparse.ArgumentParser(description='Overlap gene reporter')
parser.add_argument('gff1', type=str, help='GFF1 file')
parser.add_argument('--gff2', type=str, help='GFF2 file')
arg = parser.parse_args()

start = datetime.now()

gff = {}
with open(arg.gff1, 'r') as fp:
    for line in fp:
        line = line.strip()
        x = line.split('\t')
        idx = x[0]
        if idx not in gff:
            gff[idx] = []
        beg = x[1]
        end = x[2]
        gff[idx].append({
            'beg': int(beg),
            'end': int(end)
        })
# Read the second file if there is one
if arg.gff2:
    with open(arg.gff2, 'r') as fp:
        for line in fp:
            line = line.strip()
            x = line.split('\t')
            idx = x[0]
            if idx not in gff:
                gff[idx] = []
            beg = x[1]
            end = x[2]
            gff[idx].append({
                'beg': int(beg),
                'end': int(end)
            })
# Prepare for algorithm
events = []
for name in gff:
    for content in gff[name]:
        beg = content['beg']
        end = content['end']
        events.append({
            'position': beg,
            'priority': 0,
            'scope': (beg, end),
            'name': name,
        })
        events.append({
            'position': end,
            'priority': 1,
            'scope': (beg, end),
            'name': name,
        })
events.sort(key=lambda x: (x['position'], x['priority']))
# Algorithm
overlaps = {}
overlap = set()
for event in events:
    position = event['position']
    priority = event['priority']
    scope = event['scope']
    name = event['name']
    beg, end = scope
    index = (name, beg, end)
    if index not in overlaps:
        overlaps[index] = set()
    if priority == 0:
        # Append overlap
        for olp in overlap:
            if olp not in overlaps:
                overlaps[olp] = set()
            overlaps[olp].add((name, beg, end))
            overlaps[index].add(olp)  # Also add existing features to overlaps of current feature
        # Add overlap index
        overlap.add(index)
    elif priority == 1:
        overlap.remove(index)
# Output
output = {}
for tuple1 in overlaps:
    if tuple1 not in output:
        output[tuple1] = []
    for tuple2 in overlaps[tuple1]:
        n2, b2, e2 = tuple2
        output[tuple1].append((n2, b2, e2))
# Format the output
output_key = sorted(overlaps.keys(), key=lambda x: (x[0], x[1], x[2]))
for tuple1 in output_key:
    name1, beg1, end1 = tuple1
    if output[tuple1]:
        print(f'{name1} {beg1} {end1}')
        output[tuple1].sort(key=lambda x: (x[0], x[1], x[2]))
        for overlap_feature in output[tuple1]:
            name2, beg2, end2 = overlap_feature
            print(f'\t{name2} {beg2} {end2}')
        print()
    else:
        print(f'In {name1} {beg1} {end1} there is no overlap.\n')

end = datetime.now()
print(f"Total processing time: {end - start}")