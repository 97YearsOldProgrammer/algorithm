import argparse

parser = argparse.ArgumentParser(description='Overlap gene reporter')
parser.add_argument('gff1', type=str, help='GFF1 file')
parser.add_argument('gff2', type=str, help='GFF2 file')
arg = parser.parse_args()

def read_file(*filenames):
    events = []
    for file in filenames:
        with open(file, 'r') as fp:
            for line in fp:
                idx, beg, end = line.strip().split('\t')
                beg = int(beg)
                end = int(end)
                events.append( (idx, beg, end, 0) )
                events.append( (idx, end, beg, 1) )
    return events

events= read_file(arg.gff1, arg.gff2)
events.sort(key=lambda x: ( x[1], x[3] ))

overlaps = []
for idx, beg, end, typ in events:

    if typ == 0:
        output = (idx, beg, end)
        overlaps.append(output)

    elif typ == 1:
        output = (idx, end, beg)
        overlaps.remove(output)

        if overlaps:
            print(f'{output}--',ends="")
            print(overlaps[:])

        elif not overlaps: 
            print(f'There is no overlap for {output}')