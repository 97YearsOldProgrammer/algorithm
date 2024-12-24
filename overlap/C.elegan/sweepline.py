import argparse
import gzip
from datetime import datetime

parser = argparse.ArgumentParser(description='Overlap finder for transcripts exons/introns')
parser.add_argument('gff', type=str, help='GFF file (can be gzipped).')
arg = parser.parse_args()

start = datetime.now()

feature_types = ['exon', 'intron']
genes = {}

# Read the GFF file and parse features
with gzip.open(arg.gff, 'rt') as fp:
    for line in fp:
        line = line.strip()
        if line.startswith('#') or not line:
            continue  # Skip header and empty lines
        x = line.split('\t')
        if len(x) < 9:
            continue  # Skip invalid lines
        src = x[1]
        if src != 'WormBase':
            continue
        typ = x[2]
        if typ not in feature_types:
            continue
        chm = x[0]
        beg = int(x[3])
        end = int(x[4])
        att = x[8]
        # Extract gene/transcript name
        if ';' in att:
            att = att.split(';')[0]
        if ':' in att:
            att = att.split(':')[1]
        name = att.strip()
        if name not in genes:
            genes[name] = []
        genes[name].append({
            'type': typ,
            'beg': beg,
            'end': end
        })

# Create events for the sweep line algorithm
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
events.sort(key=lambda x: (x['position'], x['priority']))

# Initialize data structures
overlaps = {}  # Key: (name, beg, end), Value: set of overlapping features
active_features = set()

# Process events to find overlaps
for event in events:
    position = event['position']
    priority = event['priority']
    name = event['name']
    typ = event['type']
    beg, end = event['scope']
    index = (name, beg, end)
    if index not in overlaps:
        overlaps[index] = set()
    if priority == 0:
        # Feature starts; check for overlaps with active features
        for active in active_features:
            overlaps[active].add(index)
            overlaps[index].add(active)
        active_features.add(index)
    elif priority == 1:
        # Feature ends; remove from active features
        active_features.remove(index)

# Prepare and print the output
output_entries = {}
for feature_key, overlap_set in overlaps.items():
    name1, beg1, end1 = feature_key
    if feature_key not in output_entries:
        output_entries[feature_key] = set()
    for overlap_feature in overlap_set:
        name2, beg2, end2 = overlap_feature
        if (name1, beg1, end1) != (name2, beg2, end2):
            output_entries[feature_key].add((name2, beg2, end2))

# Sort the features for consistent output
sorted_features = sorted(output_entries.keys(), key=lambda x: (x[0], x[1], x[2]))

for feature in sorted_features:
    name1, beg1, end1 = feature
    overlaps_list = output_entries[feature]
    if overlaps_list:
        print(f'{name1} {beg1} {end1}')
        # Sort overlaps for consistent output
        sorted_overlaps = sorted(overlaps_list, key=lambda x: (x[0], x[1], x[2]))
        for overlap_feature in sorted_overlaps:
            name2, beg2, end2 = overlap_feature
            print(f'\t{name2} {beg2} {end2}')
        print()
    else:
        print(f'In {name1} {beg1} {end1} there is no overlap.\n')

end = datetime.now()
print(f"Total processing time: {end - start}")

