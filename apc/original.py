import itertools

def short_intron(dons, accs, min):
	for d, a in zip(dons, accs):
		intron_length = a - d + 1
		if intron_length < min: return True
	return False

def short_exon(dons, accs, seqlen, flank, min):

	# first exon
	exon_beg = flank + 1
	exon_end = dons[0] -1
	exon_len = exon_end - exon_beg + 1
	if exon_len < min: return True

	# last exon
	exon_beg = accs[-1] + 1
	exon_end = seqlen - flank + 1
	exon_len = exon_end - exon_beg + 1
	if exon_len < min: return True

	# interior exons
	for i in range(1, len(dons)):
		exon_beg = accs[i-1] + 1
		exon_end = dons[i] - 1
		exon_len = exon_end - exon_beg
		if exon_len < min: return True
	return False

def gtag_sites(seq, flank, minex):
	dons = []
	accs = []
	for i in range(flank + minex, len(seq) -flank -minex -1):
		if seq[i:i+2]   == 'GT': dons.append(i)
		if seq[i-1:i+1] == 'AG': accs.append(i)
	return dons, accs

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
			beg = int(f[3]) -1
			end = int(f[4]) -1
			if gtag:
				if seq[beg:beg+2]   != 'GT': continue
				if seq[end-1:end+1] != 'AG': continue
			dond[beg] = True
			accd[end] = True

	dons = list(dond.keys())
	accs = list(accd.keys())
	dons.sort()
	accs.sort()

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
	tx['exons'].append((beg, dons[0] -1))
	for i in range(1, len(dons)):
		a = accs[i-1] + 1
		b = dons[i] -1
		tx['exons'].append((a, b))
	tx['exons'].append((accs[-1] +1, end))

	return tx

def all_possible(seq, minin, minex, maxs, flank, gff=None):

	if gff: dons, accs = gff_sites(seq, gff)
	else:   dons, accs = gtag_sites(seq, flank, minex)

	info = {
		'trials' : 0,
		'donors': len(dons),
		'acceptors': len(accs),
		'short_intron': 0,
		'short_exon': 0,
	}

	isoforms = []
	sites = min(len(dons), len(accs), maxs)
	for n in range(1, sites+1):
		for dsites in itertools.combinations(dons, n):
			for asites in itertools.combinations(accs, n):
				info['trials'] += 1

				# sanity checks
				if short_intron(dsites, asites, minin):
					info['short_intron'] += 1
					continue

				if short_exon(dsites, asites, len(seq), flank, minex):
					info['short_exon'] += 1
					continue

				# create isoform and save
				tx = build_mRNA(seq, flank, len(seq) -flank -1, dsites, asites)
				isoforms.append(tx)

	return isoforms, info


"""                     test area                    """


seq = '''AGGTTTTCAGTCTGCGTATCTAAATCAGCATTAATAACAATGAGGCGAGGATAACAGCAA
TCGGGCCTCTAGTCCTTTTTGCTGGCCAAAGTGAACGTCAGCTAGGTAACGTTGCGCATT
CAAAGACGAGTAGTGTCGCGAAACAAATAGTTATCAGTGCTACAATTACATGCGCAACAA
TTGAGGAGTTAAGGCACTACCACGGAACAGAGACCGTGTGCCAGGAGAAGTGGTACGATT
GCTCCTGACAGTTGCCGACTCCGGGCGCATTGACTGTCACTTCGAAGTCATGGCA
'''

isoforms, info = all_possible(seq, 20, 20, 5, 0)
print(len(isoforms), info)