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

        # the exit mechanism
        if i == 2:
            # check last exon
            if len(seq) - flank - sol[-1] + 1 < minex: 
                info['short_exon'] += 1
                return
            intron1.append(sol[:])
            return

        if i == 0:
        # we pick the donor first as parental node
            for n1 in dons:
                info['trails'] += 1
                
                # sanity check
                if not sol and n1 - 1 - flank < minex:
                    info['short_exon'] += 1
                    continue
                elif sol and n1 < sol[-1] + minex + 2:
                    info['short_exon'] += 1 
                    continue

                # algorithm
                sol.append(n1)
                backtrack(i + 1)
                sol.pop()

        # where we pick the acceptor as parental node
        else:
            for n2 in accs:
                info['trails'] += 1
                # sanity check
                if sol and n2 < sol[-1]: continue
                if sol and n2 < sol[-1] + minin - 1: 
                    info['short_intron'] += 1
                    continue
                # algorithm
                sol.append(n2)
                backtrack(i + 1)
                sol.pop()

    backtrack(0)

    sol = []
    isoforms = []

    def dynamic_programming(n):

        if n == maxs: 
            tx = build_mRNA(seq, flank, len(seq) - flank - 1, sol[:] )
            isoforms.append(tx)
            return
        
        if sol:
            tx = build_mRNA(seq, flank, len(seq) - flank - 1, sol[:] )
            isoforms.append(tx)
        
        for i1 in intron1:
            info['trails'] += 1
            if sol and i1[0] <= sol[-1]: continue
            if sol and i1[0] < sol[-1] + minex + 2:
                info['short_exon'] += 1
                continue
            sol.extend(i1)
            dynamic_programming(n+1)
            sol.pop()
            sol.pop()

    dynamic_programming(0)
    return isoforms, info