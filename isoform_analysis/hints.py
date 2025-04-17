import subprocess

def run(model, fasta):
    '''
    i know it's stupid, hardcode path, maybe better plan later on
    '''

    model_files = [
        "../isoform_analysis/models/don.pwm",
        "../isoform_analysis/models/acc.pwm",
        "../isoform_analysis/models/exon.mm",
        "../isoform_analysis/models/intron.mm",
        "../isoform_analysis/models/exon.len",
        "../isoform_analysis/models/intron.len"
    ]
    
    cmd = [model, fasta] + model_files
    result = subprocess.run(cmd, check=True, text=True, capture_output=True)
    
    return result.stdout

def parse(output):
    '''
    Parse EDHMM output to extract donor and acceptor site probabilities
    '''
    
    dons = []
    accs = []
    switch = 0
    
    lines = output.strip().split('\n')

    for line in lines:
        # this switch logic is odd, lit can turn into anything better
        if      line == "DONS": 
            switch = 0
            continue

        elif    line == "ACCS":
            switch = 1
            continue
        
        line = line.strip().split('\t')
        pos = int  (line[0])
        val = float(line[1])

        if   switch == 0: dons.append( (pos, val) )
        elif switch == 1: accs.append( (pos, val) )
    
    return dons, accs

def gapstats(sites):
    '''
    Gap statistic to find significant splice sites
    '''
    sites = sorted(sites, key=lambda x: x[1], reverse=True)
    max = 0
    idx = 0
    
    for i in range( len(sites) - 1 ):
        val  = sites[i][1]
        val2 = sites[i+1][1]
        gap  = val / val2

        if gap > max:
            max = gap
            idx = i
    
    return sorted( [site[0] for site in sites[:idx+1]] )

def smoothed_gapstats(sites, k:int = 5):
    '''
    smoothed gap statistics 
    '''
    sites= sorted(sites, key=lambda x: x[1], reverse=True)
    max = 0
    idx = 0

    for i in range( len(sites) - k ):
        avg  = sum(site[1] for site in sites[i:i+k])     / k
        avg2 = sum(site[1] for site in sites[i+1:i+k+1]) / k
        gap  = avg/avg2 

        if gap > max:
            max = gap
            idx = i + k // 2
    
    return sorted( [site[0] for site in sites[:idx+1]] )

def percentile(sites, percentile:int = 25):
    '''
    filter splice sites by percentile
    '''

    sites  = sorted(sites, key=lambda x:[1], reverse=True)
    cutoff = sites[int ( len(sites) * (1 - percentile/100) )][1]
    sites  = [site for site in sites if site[1] >= cutoff]
    return sorted([site[0] for site in sites])