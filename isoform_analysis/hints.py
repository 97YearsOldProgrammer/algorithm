import os
import subprocess
import argparse
import korflab

def run(model, fasta):

    '''
    Run EDHMM model on a single sequence file and return the output
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
    swtich = 0
    
    for line in output.strip().split('\n'):
        
        if      line == "DONS": 
            switch = 0
            continue
        
        elif    line == "ACCS":
            swtich = 1
            continue
        
        if   swtich == 0: dons.append[float(line)]
        elif swtich == 1: accs.append[float(line)]
    
    return dons, accs

def analyze_splice_sites(dons, accs):
    """
    Analyze splice sites using gap statistic to identify significant sites
    
    Args:
        dons: Dictionary of donor site probabilities
        accs: Dictionary of acceptor site probabilities
    
    Returns:
        Tuple of (significant_dons, significant_accs) dictionaries
    """
    # Find significant donor sites using gap statistic
    significant_dons = find_significant_sites(dons)
    
    # Find significant acceptor sites using gap statistic
    significant_accs = find_significant_sites(accs)
    
    return significant_dons, significant_accs

def gstats(sites):
    '''
    Implement gap statistic to find significant splice sites
    '''
    sorted_sites = sorted(sites.items(), key=lambda x: x[1], reverse=True)
    
    # If we have fewer than 2 sites, return all of them
    if len(sorted_sites) < 2:
        return dict(sorted_sites)
    
    max_gap_ratio = 0
    cutoff_idx = 0
    
    # Find the largest gap between consecutive values
    for i in range(len(sorted_sites) - 1):
        if sorted_sites[i][1] <= 0:
            continue
            
        gap_ratio = sorted_sites[i][1] / sorted_sites[i+1][1]
        if gap_ratio > max_gap_ratio:
            max_gap_ratio = gap_ratio
            cutoff_idx = i
    
    # Return only sites before the significant gap
    return dict(sorted_sites[:cutoff_idx + 1])

def gtag_sites(seq, flank, minex):
    
    dons = []
    accs = []
    for i in range(flank + minex, len(seq) - flank - minex - 1):
        if seq[i:i+2] == 'GT': dons.append(i)
        if seq[i-1:i+1] == 'AG': accs.append(i)
    return dons, accs



if __name__ == "__main__":
    main()