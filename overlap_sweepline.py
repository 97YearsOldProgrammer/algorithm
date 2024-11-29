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

# List to hold all events
events = []
# Dictionary to hold feature details
feature_details = {}

# Build events and feature_details
feature_counter = 0  # Unique ID for each feature
for trxn in genes:
    for feature in genes[trxn]:
        feature_id = feature_counter
        feature_counter += 1
        feature['name'] = trxn
        feature_details[feature_id] = feature
        events.append({
            'position': feature['beg'],
            'event_type': 0,  # 0 for start
            'feature_id': feature_id
        })
        events.append({
            'position': feature['end'],
            'event_type': 1,  # 1 for end
            'feature_id': feature_id
        })

# Sort the events
events.sort(key=lambda x: (x['position'], x['event_type']))

# Initialize variables for overlap detection
active_features = {}  # key: feature_id, value: feature_details
overlaps = []

# Process events to find overlaps
for event in events:
    position = event['position']
    feature_id = event['feature_id']
    event_type = event['event_type']

    current_feature = feature_details[feature_id]

    if event_type == 0:  # Start of a feature
        # Before adding the current feature, check for overlaps with active features
        for other_feature_id, other_feature in active_features.items():
            # Check if current_feature overlaps with other_feature
            overlap_start = max(current_feature['beg'], other_feature['beg'])
            overlap_end = min(current_feature['end'], other_feature['end'])
            if overlap_start < overlap_end:
                # There is an overlap
                overlaps.append({
                    'overlap_start': overlap_start,
                    'overlap_end': overlap_end,
                    'features': [current_feature, other_feature]
                })
        # Add the current feature to active features
        active_features[feature_id] = current_feature
    else:  # End of a feature
        # Remove the feature from active features
        active_features.pop(feature_id, None)

# Output the overlaps with feature details
for overlap in overlaps:
    overlap_start = overlap['overlap_start']
    overlap_end = overlap['overlap_end']
    f1, f2 = overlap['features']
    print(f"{f1['name']} {f1['type']}  {f1['beg']} {f1['end']} overlaps with "
          f"{f2['name']} {f2['type']}  {f2['beg']} {f2['end']} ")


end = datetime.now()
print(f"Total processing time: {end - start}")
