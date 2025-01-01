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

def build_mRNA(seq, beg, end, dons, accs):
	assert(beg <= end)
	assert(len(dons) == len(accs))

	tx = {
		'seq': seq,
		'beg': beg,
		'end': end,
		'exons': [],
		'introns': [],
		'score': 0
	}

	if len(dons) == 0:
		tx['exons'].append((beg, end))
		return tx

	# introns
	for a, b in zip(dons, accs):
		tx['introns'].append((a, b))

	# exons
	tx['exons'].append((beg, dons[0] - 1))
	for i in range(1, len(dons)):
		a = accs[i-1] + 1
		b = dons[i]   - 1
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

    don = []
    acc = []
    isoforms = []

    def backtrack(i):

        if i == maxs * 2: return

        if i % 2 == 0:

            for ds in dons:
                info['trails'] += 1

                # fast check
                if don and ds <= acc[-1]: continue
                # sanity check
                else: 
                    # check first exon
                    if ds - 1 - flank < minex:
                        info['short_exon'] += 1
                        continue

                # check short_exon
                if acc and ds < acc[-1] + minex + 2:
                    info['short_exon'] += 1
                    continue

                # algorithm
                don.append(ds)
                backtrack(i + 1)
                don.pop()

        else:

            for ac in accs:
                info['trails'] += 1

                # fast check
                if ac <= don[-1]: continue
                # check short_intron
                if ac < don[-1] + minin - 1:
                    info['short_intron'] += 1
                    continue

                acc.append(ac)
                # check last exon
                if len(seq) - flank - ac + 1 >= minex:
                    # creating isoform and store
                    tx = build_mRNA(seq, flank, len(seq) - flank - 1, don, acc )
                    isoforms.append(tx)
                    backtrack(i + 1)
                    acc.pop()
                else: 
                    info['short_exon'] += 1
                    acc.pop()
                    continue
    
    backtrack(0)
    return isoforms, info