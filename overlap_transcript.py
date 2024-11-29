import argparse
import optool
import gzip
from datetime import datetime
from collections import defaultdict

start = datetime.now()

parser = argparse.ArgumentParser(description='Overlap gene reporter')
parser.add_argument('gff', type=str, help='GFF file (gzip file)')
arg = parser.parse_args()

feature_types = ['exon', 'intron', 'five_prime_UTR', 'three_prime_UTR']
transcript_types = [
    'mRNA', 'piRNA', 'snoRNA', 'miRNA', 'pre_miRNA',
    'pseudogenic_transcript', 'tRNA', 'pseudogenic_tRNA',
    'ncRNA', 'circular_ncRNA', 'nc_primary_transcript',
    'transcript']

genes = defaultdict(                            #chromosome
    lambda: defaultdict(                        #gene
        lambda: {       
            'type': '',                         #gene biotype
            'strand': '',                       #gene on which strand
            'transcripts': defaultdict(         #different transcript
                lambda: {                       #transcript
                    'type': '',                 #transcript biotype
                    'features': [],             #contain multiple introns,exons,beg, end
                    'feature_set': set()        #remove repetitive
                }
            )
        }
    )
)

with gzip.open(arg.gff, 'rt') as fp:
    current_gene = None
    current_chm = None

    for line in fp:
        x = line.strip().split('\t')
        src = x[1]
        if src != 'WormBase': continue
        chm = x[0]
        typ = x[2]
        beg = x[3]
        end = x[4]
        pol = x[6]
        att = x[8]

        if typ == 'gene':
            current_gene = optool.gene_namer(att, 'sequence_name')
            gene_biotype = optool.gene_namer(att, 'biotype')
            genes[chm][current_gene]['strand'] = pol
            genes[chm][current_gene]['type'] = gene_biotype
            current_chm = chm # for chm
            continue

        if current_gene:
            # transcript type
            if 'ID=' in att:
                transcript = optool.transcript_namer(att, 'ID')
                if transcript:
                    transcript_biotype = optool.gene_namer(att, 'biotype')
                    # if there is nothing we could find for biotype
                    if not transcript_biotype:
                        transcript_biotype = typ  
                    genes[current_chm][current_gene]['transcripts'][transcript]['type'] = transcript_biotype
            # if it is a intron, exons
            if typ in feature_types:
                feature_info = {
                    'type': typ,
                    'start': int(beg),
                    'end': int(end),
                    'attributes': att
                }

                def get_parents(attributes, key):
                    attributes = attributes.strip().split(';')
                    for attribute in attributes:
                        attribute = attribute.strip()
                        if attribute.startswith(f'{key}='):
                            names = attribute.split('=', 1)[1]
                            names = names.strip('"').split(',')
                            parents = []
                            prefixes = ['Transcript:', 'Gene:', 'Pseudogene:', 'Sequence:']
                            for name in names:
                                name = name.strip()
                                for prefix in prefixes:
                                    if name.startswith(prefix):
                                        name = name[len(prefix):]
                                        break
                                parents.append(name)
                            return parents
                    return []

                parents = get_parents(att, 'Parent')
                if not parents:
                    parents = [current_gene]

                feature_key = (typ, int(beg), int(end))
                gene_entry = genes[current_chm][current_gene]
                transcripts = gene_entry['transcripts']

                for parent in parents:
                    parent = parent.strip()
                    if parent in transcripts:
                        transcript_entry = transcripts[parent]
                    else:
                        # If parent transcript not found, create it
                        transcript_entry = transcripts[parent]
                        # Optionally set the transcript type to the gene's type
                        transcript_entry['type'] = gene_entry['type']

                    if feature_key not in transcript_entry['feature_set']:
                        transcript_entry['features'].append(feature_info)
                        transcript_entry['feature_set'].add(feature_key)

for chm in genes:
    print(f"Chromosome: {chm}")
    for gene_name in genes[chm]:
        gene_entry = genes[chm][gene_name]
        strand = gene_entry['strand']
        gene_type = gene_entry['type']
        print(f"  Gene: {gene_name} (Type: {gene_type}, Strand: {strand})")
        transcripts = gene_entry['transcripts']
        for transcript_name in transcripts:
            transcript_entry = transcripts[transcript_name]
            transcript_type = transcript_entry['type']
            print(f"    Transcript: {transcript_name} (Type: {transcript_type})")
            # sort by the beg position
            sorted_features = sorted(transcript_entry['features'], key=lambda f: f['start'])
            for feature in sorted_features:
                print(f"      Feature: {feature['type']} ({feature['start']}-{feature['end']})")

end = datetime.now()
print(end-start)