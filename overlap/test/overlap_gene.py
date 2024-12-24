import argparse
import optool
from datetime import datetime

start = datetime.now()

parser = argparse.ArgumentParser(description='overlap gene reporter')
parser.add_argument('gff', type=str, help='GFF file (gzip file)')
arg = parser.parse_args()
          
gene = optool.overlap_gene_gff(arg.gff)

header = True
for chm in gene:

    gene_list = gene[chm]
    gene_list = sorted(gene_list, key=lambda x: int(x[0])) 
    
    gene_indices = {}
    idx = 1
    for beg, end, pol, seq in gene_list:
        if seq not in gene_indices:
            gene_indices[seq] = idx
            idx += 1

    events = []
    for beg, end, pol, seq in gene_list:
        idx = gene_indices[seq]
        events.append({
            'position': beg,
            'polarity': pol,
            'index': idx,
            'type': 'beg',
            'name': seq
        })
        events.append({
            'position': end,
            'polarity': pol,
            'index': idx,
            'type': 'end',
            'name': seq
        })

    events = sorted(events, key=lambda x: int(x['position']))

    genes = set() 
    overlaps = []

    for event in events:
        typ = event['type']
        idx = event['index']
        pol = event['polarity']
        name = event['name']

        if typ == 'beg':
            for active_idx, active_pol, active_name in genes:
                # Define overlap type based on polarities
                if pol == active_pol:
                    overlap_typ = 'Same-Strand'
                else:
                    overlap_typ = 'Opposite-Strand'
                overlaps.append((idx, name, pol, active_idx, active_name, active_pol, overlap_typ))
            genes.add((idx, pol, name))
        
        elif typ == 'end':
            genes.remove((idx, pol, name))
    
    overlaps = sorted(overlaps, key=lambda x: (x[0], x[3]))
    if header:
        header_text = f'{"Chromosome":<10}\t{"Gene1":<15}\t{"Pol1":<5}\t{"Idx1":<5}\t{"Gene2":<15}\t{"Pol2":<5}\t{"Idx2":<5}\t{"Overlap_Type"}'
        print(header_text)
        header = False

    for idx1, name1, pol1, idx2, name2, pol2, overlap_typ in overlaps:
        print(f'{chm:<10}\t{name1:<15}\t{pol1:<5}\t{idx1:<5}\t{name2:<15}\t{pol2:<5}\t{idx2:<5}\t{overlap_typ}')

end = datetime.now()
print(f"Total processing time: {end - start}")
# print(f'Time taken: {end - start}')

