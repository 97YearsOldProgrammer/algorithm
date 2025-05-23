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

def all_possible(seq, minin, minex, maxs, flank, gff=None):
    
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
    intron1  = []

    def backtrack(i):

        if i == 2: return
 
        if i % 2 == 0:

            for ds in dons:
                info['trails'] += 1

                # check first exon
                if not sol and ds - flank < minex:
                    info['short_exon'] += 1
                    continue
                # check interior exon
                if sol and ds < sol[-1] + minex + 1:
                    info['short_exon'] += 1 
                    continue

                # algorithm
                sol.append(ds)
                backtrack(i + 1)
                sol.pop()

        else:
            
            for ac in accs:
                info['trails'] += 1

                # fast check
                if ac <= sol[-1]: continue
                # sanity check
                if ac < sol[-1] + minin - 1: 
                    info['short_intron'] += 1
                    continue

                sol.append(ac)
                # check last exon
                if len(seq) - flank - ac - 1 >= minex:
                    intron1.append( sol[:] )
                    backtrack(i + 1)
                    sol.pop()
                else: 
                    info['short_exon'] += 1
                    sol.pop()
                    continue

    backtrack(0)

    sol = []
    isoforms = []

    def dynamic_programming(n):

        if n == maxs: return

        for i1 in intron1:
            info['trails'] += 1
            # sanity check
            if sol and i1[0] <= sol[-1]: continue
            if sol and i1[0] < sol[-1] + minex + 1:
                info['short_exon'] += 1
                continue

            sol.extend(i1)
            # create isoforms and store
            tx = build_mRNA(seq, flank, len(seq) - flank - 1, sol[:] )
            isoforms.append(tx)
            # backtrack
            dynamic_programming(n+1)
            sol.pop()
            sol.pop()

    dynamic_programming(0)
    return isoforms, info