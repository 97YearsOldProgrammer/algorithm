import argparse
import gzip
import optool

parser = argparse.ArgumentParser(description='Overlap finder for two genes')
parser.add_argument('gff', type=str, help='GFF file (can be gzipped).')
parser.add_argument('gene1', type=str, help='Gene1 name')
parser.add_argument('gene2', type=str, help='Gene2 name')
arg = parser.parse_args()

feature_types = ['exon', 'intron']
gene1_active = False
gene2_active = False

genes = {arg.gene1: [], arg.gene2: []}

with gzip.open(arg.gff, 'rt') as fp:
    for line in fp:
        x = line.strip().split('\t')
        src = x[1]
        if src != 'WormBase':
            continue
        typ = x[2]
        beg = int(x[3])
        end = int(x[4])
        att = x[8]

        if typ == 'gene':
            gene_name = optool.gene_namer(att, 'sequence_name')
            if gene_name == arg.gene1:
                gene1_active = True
                gene2_active = False
            elif gene_name == arg.gene2:
                gene2_active = True
                gene1_active = False
            else:
                gene1_active = False
                gene2_active = False

        if gene1_active or gene2_active:
            if typ in feature_types:
                if 'Parent=' in att:
                    name = optool.gene_namer(att, 'Parent')
                    if name.startswith("Transcript:"):
                        name = name.replace("Transcript:", "", 1)
                    gene_key = arg.gene1 if gene1_active else arg.gene2
                    genes[gene_key].append((beg, end, typ, name, gene_key))

# Now, process the intervals to find overlaps
features_gene1 = sorted(genes[arg.gene1], key=lambda x: x[0])
features_gene2 = sorted(genes[arg.gene2], key=lambda x: x[0])

overlaps = []

j = 0
for feature1 in features_gene1:
    beg1, end1, typ1, trans1, gene1_name = feature1

    # Skip features in gene2 that end before feature1 starts
    while j < len(features_gene2) and features_gene2[j][1] < beg1:
        j += 1

    k = j
    while k < len(features_gene2) and features_gene2[k][0] <= end1:
        feature2 = features_gene2[k]
        beg2, end2, typ2, trans2, gene2_name = feature2

        # Check for overlap
        overlap_beg = max(beg1, beg2)
        overlap_end = min(end1, end2)
        if overlap_beg <= overlap_end:
            overlaps.append({
                'overlap_start': overlap_beg,
                'overlap_end': overlap_end,
                'gene1': gene1_name,
                'transcript1': trans1,
                'feature_type1': typ1,
                'gene2': gene2_name,
                'transcript2': trans2,
                'feature_type2': typ2
            })
        k += 1

if overlaps:
    header = (f"{'Gene1':<15}{'Transcript1':<20}{'Feature1':<10}"
                   f"{'Gene2':<15}{'Transcript2':<20}{'Feature2':<10}{'Overlap_Location'}")
    print(header)

    for overlap in overlaps:
        location = f"{overlap['overlap_start']}-{overlap['overlap_end']}"
        print(f"{overlap['gene1']:<15}{overlap['transcript1']:<20}{overlap['feature_type1']:<10}"
            f"{overlap['gene2']:<15}{overlap['transcript2']:<20}{overlap['feature_type2']:<10}{location}")
else: print('Sadly, there is no overlap between this two gene.')