import argparse
import gzip
import optool
from datetime import datetime

# Argument parsing
parser = argparse.ArgumentParser(description='Overlap finder for transcripts exons/introns')
parser.add_argument('gff', type=str, help='GFF file (can be gzipped).')
arg = parser.parse_args()

start = datetime.now()

# Specify feature types of interest
feature_types = ['exon', 'intron']

# Dictionary to store transcripts and their features
genes = {}

# Read the GFF file and build the genes dictionary
with gzip.open(arg.gff, 'rt') as fp:
    for line in fp:
        x = line.strip().split('\t')
        if len(x) < 9:
            continue  # Skip malformed lines
        src = x[1]
        if src != 'WormBase': 
            continue
        typ = x[2]
        if typ not in feature_types: 
            continue
        chm = x[0]
        beg = x[3]
        end = x[4]
        att = x[8]
        if 'Parent=' in att:
            name = optool.gene_namer(att, 'Parent')
            if name.startswith("Transcript:"):
                name = name.replace("Transcript:", "", 1)
                if name not in genes:
                    genes[name] = []
                genes[name].append({
                    'type': typ,
                    'beg': int(beg),
                    'end': int(end)
                })

events = []
feature_details = {}

idx = -1  
for trxn in genes:
    for feature in genes[trxn]:
        idx += 1
        feature['name'] = trxn
        feature_details[id] = feature
        events.append({
            'position': feature['beg'],
            'event_type': 0, 
            'id': id
        })
        events.append({
            'position': feature['end'],
            'event_type': 1, 
            'id': id
        })

events.sort(key=lambda x: (x['position'], x['event_type']))

active_features = {}  # key: feature_id, value: feature_details
overlaps = []

for event in events:
    id = event['id']
    position = event['position']
    event_type = event['event_type']

    current_feature = feature_details[id]

    if event_type == 0:  
        for other_feature_id, other_feature in active_features.items():
            overlap_start = max(current_feature['beg'], other_feature['beg'])
            overlap_end = min(current_feature['end'], other_feature['end'])
            if overlap_start < overlap_end:
                overlaps.append({
                    'overlap_start': overlap_start,
                    'overlap_end': overlap_end,
                    'features': [current_feature, other_feature]
                })
        active_features[id] = current_feature
    else:  
        active_features.pop(id, None)

for overlap in overlaps:
    overlap_start = overlap['overlap_start']
    overlap_end = overlap['overlap_end']
    f1, f2 = overlap['features']
    print(f"{f1['name']} {f1['type']}  {f1['beg']} {f1['end']} overlaps with "
          f"{f2['name']} {f2['type']}  {f2['beg']} {f2['end']} ")


end = datetime.now()
print(f"Total processing time: {end - start}")
