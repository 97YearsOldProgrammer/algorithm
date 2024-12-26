def gff_sites(seq, gff, gtag=True):
    dond = {}
    accd = {}
    with open(gff) as fp:
        for line in fp:
            if line.startswith('#'): continue
            f = line.split()
            if len(f) < 8: continue
            if f[2] != 'intron': continue
            if f[6] != '+': continue
            beg = int(f[3]) - 1
            end = int(f[4]) - 1
            if gtag:
                if seq[beg:beg+2] != 'GT': continue
                if seq[end-1:end+1] != 'AG': continue
            dond[beg] = True
            accd[end] = True

    dons = list(dond.keys())
    accs = list(accd.keys())
    dons.sort()
    accs.sort()

    return dons, accs

def gtag_sites(seq, flank, minex):
    dons = []
    accs = []
    for i in range(flank + minex, len(seq) - flank - minex - 1):
        if seq[i:i+2] == 'GT': dons.append(i)
        if seq[i-1:i+1] == 'AG': accs.append(i)
    return dons, accs

def build_mRNA(seq, beg, end, res):
    assert(beg <= end)
    assert(len(res) % 2 == 0)

    tx = {
        'seq': seq,
        'beg': beg,
        'end': end,
        'exons': [],
        'introns': [],
        'score': 0
    }

    dons = []
    accs = []

    for i, index in enumerate(res):
        if i % 2 == 0: dons.append(index)
        else: accs.append(index)

    # introns
    for a, b in zip(dons, accs):
        tx['introns'].append((a, b))

    # exons
    tx['exons'].append((beg, dons[0] - 1))
    for i in range(1, len(dons)):
        a = accs[i-1] + 1
        b = dons[i] - 1
        tx['exons'].append((a, b))
    tx['exons'].append((accs[-1] + 1, end))

    return tx

def single_possible(seq, minin, minex, maxs, flank, gff=None):
    
    if gff: dons, accs = gff_sites(seq, gff)
    else:   dons, accs = gtag_sites(seq, flank, minex)

    info = {
            'trails': 0,
            'donors': len(dons),
            'acceptors': len(accs),
            'short_intron': 0,
            'short_exon': 0,
    }

    sol = []
    introns  = []
    isoforms = []

    def backtrack(i):

        if i == maxs*2:
            # check last exon
            if len(seq) - flank - sol[-1] + 1 < minex: 
                info['short_exon'] += 1
                return
            # create isoform and save
            introns.append(sol[:])
            tx = build_mRNA(seq, flank, len(seq) - flank - 1, sol[:] )
            isoforms.append(tx)
            return

        if (maxs*2 - i) % 2 == 0:
            for n1 in range(len(dons)):
                info['trails'] += 1

                # sanity check
                if not sol:
                    if dons[n1] - 1 - flank < minex:
                        info['short_exon'] += 1
                        continue
                elif sol and dons[n1] < sol[-1] + minex + 2:
                    info['short_exon'] += 1 
                    continue

                # algorithm
                sol.append(dons[n1])
                backtrack(i + 1)
                sol.pop()

        else:
            for n2 in range(len(accs)):
                info['trails'] += 1

                # sanity check
                if sol and accs[n2] < sol[-1]: continue
                if sol and accs[n2] < sol[-1] + minin - 1: 
                    info['short_intron'] += 1
                    continue
                # algorithm
                sol.append(accs[n2])
                backtrack(i + 1)
                sol.pop()

    backtrack(0)
    return isoforms, info, introns

def all_possible(seq, minin, minex, maxs, flank, gff=None):

    isoforms = []
    info = {
        'trails': 0,
        'donors': 0,
        'acceptors': 0,
        'short_intron': 0,
        'short_exon': 0,
    }
    introns = dict()

    def dynamic_programming(seq, minin, minex, i, flank):
            introns[i] = []

            for inner_intron in introns[i - 2]:
                info['trails'] += 1
                # innner_intron check
                if inner_intron[0] - flank - 3 < minex * 2 + minin: 
                    info['short_intron'] += 1
                    info['short_exon']   += 2
                    continue
                if len(seq) - flank - inner_intron[-1] - 3 < minex * 2 + minin:
                    info['short_intron'] += 1
                    info['short_exon']   += 2
                    continue

                for outer_intron in introns[2]:
                    info['trails'] += 1
                    # sanity check
                    if inner_intron[0]  <= outer_intron[1] : continue
                    if inner_intron[-1] >= outer_intron[2] : continue
                    # outer_intron check
                    if inner_intron[0] - outer_intron[1] - 2 < minex: 
                        info['short_exon'] += 1
                        continue
                    if outer_intron[2] - inner_intron[-1] - 2 < minex: 
                        info['short_exon'] += 1
                        continue

                    # output
                    output = outer_intron[0:2] + inner_intron[:] + outer_intron[2:]

                    # collect data for further dynamic programming
                    introns[i].append(output)

                    # create isoform and save
                    tx     = build_mRNA(seq, flank, len(seq) - flank - 1, output )
                    isoforms.append(tx)

    for i in range(1, maxs+1):

        if i == 1 or i == 2: 
            isoform, infos, intron = single_possible(seq, minin, minex, i, flank, gff)
            # collect intron
            introns[i] = []
            introns[i].extend(intron)
            # collect information
            isoforms = isoforms + isoform
            info['trails']       += infos['trails']
            info['short_intron'] += infos['short_intron']
            info['short_exon']   += infos['short_exon']
            info['donors']        = infos['donors']
            info['acceptors']     = infos['acceptors']

        else: dynamic_programming(seq, minin, minex, i, flank)
        
    return isoforms, info


"""                     test area                    """


seq = '''AGGTTTTCAGTCTGCGTATCTAAATCAGCATTAATAACAATGAGGCGAGGA
TCGGGCCTCTAGTCCTTTTTGCTGGCCAAAGTGAACGTCAGCTAGGTAACGTTGCGCATT
CAAAGACGAGTAGTGTCGCGAAACAAATAGTTATCAGTGCTACAATTACATGCGCAACAA
TTGAGGAGTTAAGGCACTACCACGGAACAGAGACCGTGTGCCAGGAGAAGTGGTACGATT
GCTCCTGACAGTTGCCGACTCCGGGCGCATTGACTGTCACTTCGAAGTCATGGCAATCAG
TCGGGCCTCTAGTCCTTTTTGCTGGCCAAAGTGAACGTCAGCTAGGTAACGTTGCGCATT
CAAAGACGAGTAGTGTCGCGAAACAAATAGTTATCAGTGCTACAATTACATGCGCAACAA
TTGAGGAGTTAAGGCACTACCACGGAACAGAGACCGTGTGCCAGGAGAAGTGGTACGATT
GCTCCTGACAGTTGCCGACTCCGGGCGCATTGACTGTCAACGTCAGCTAGGTAACGTTGC
CAAAGACGAGTAGTGTCGCGAAACAAATAGTTATCAGTGCTACAATTACATGCGCAACAA
TTGAGGAGTTAAGGCACTACCACGGAACAGAGACCAGAGACCGTGTGCCAGGAGAAGTGG
GCTCCTGACAGTTGCCGACTCCGGGCGCATTGACTGTCACTTCGAAGTCATGGCAATCAG
TCGGGCCTCTAGTCCTTTTTGCTGGCCAAAGTGAACGTCAGCTAGGTAACGTTGCGCATT
CAAAGACGAGTAGTGTCGCGAAACAAATAGTTATCAGTGCTACAATTACATGCGCAACAA
TTGAGGAGTTAAGGCACTACCACGGAACAGAGACCGTGTGCCAGGAGAAGTGGTACGATT
GCTCCTGACAGTTGCCGACTCCGGGCGCATTTTGACTGTCAACGTCAGCTAGGTAACGTT
CAAAGACGAGTAGTGTCGCGAAACAAATAGTTATCAGTGCTACAATTACATGCGCAACAA
TTGAGGAGTTAAGGCACTACCACGGAACAGAGACCAGAGACCGTGTGCCAGGAGAAGTGG
GCTCCTGACAGTTGCCGACTCCGGGCGCATTGACTGTCACTTCGAAGTCATGGCAATCAG
TCGGGCCTCTAGTCCTTTTTGCTGGCCAAAGTGAACGTCAGCTAGGTAACGTTGCGCATT
CAAAGACGAGTAGTGTCGCGAAACAAATAGTTATCAGTGCTACAATTACATGCGCAACAA
TTGAGGAGTTAAGGCACTACCACGGAACAGTCGCGAAACAAATAGTTATCAGTGCTACAA
TTGAGGAGTTAAGGCACTACCACGGAACAGAGACCGTGTGCCAGGAGAAGTGGTACGATT
GCTCCTGACAGTTGCCGACTCCGGGCGCATTTTGACTGTCAACGTCAGCTAGGTAACGTT
CAAAGACGAGTAGTGTCGCGAAACAAATAGTTATCAGTGCTACA'''

isoforms, info = all_possible(seq, 25, 35, 3, 100)
print( len(isoforms), info )