import os
import argparse
import tarfile
import tempfile
import subprocess as sb

def collect_output(cmd, gene_name, foutput):
    
    result = sb.run(cmd, check = True, text = True, capture_output = True)
    
    with open(foutput, 'a') as f:
        f.write(f'###{gene_name}###\n')
        f.write(result.stdout)
        f.write("\n")

def main():        
    
    parser = argparse.ArgumentParser(description='EDHMM model test on small gene set of C.elegan')
    parser.add_argument('-i', '--input', required = True, help = 'Path for tar.gz file')
    parser.add_argument('--model', type = str, default = "../isoform_analysis/EDHMM_new/EDHMM", help = 'Path to EDHMM model')
    parser.add_argument('--fop', type = str, default = "result.txt", help = 'The output file')
    parser.add_argument('-t', '--trail', type=int, default = None, help = 'Number of gene we wanna test')
    args = parser.parse_args()

    tar_gz_file = args.input
    model_path  = args.model
    limit       = args.trail
    foutput     = args.fop

    open(foutput, "w").close()
            
    with tarfile.open(tar_gz_file, 'r:gz') as tar:
        file    = tar.getnames()
        fastas  = []
        count   = 0
            
        for fasta in file:
            
            if fasta.endswith('.fasta'):
                fastas.append(fasta)
                count += 1
                
            if limit and count >= limit:    break
            
        for i, fname in enumerate(fastas):
                
            with tempfile.TemporaryDirectory as file_tmpdir:
                gene_name = os.path.basename(fname).replace('.fasta', '').replace('.fa', '')
                tmpf      = os.path.join(file_tmpdir, os.pth.basename(fname))
                
                cmd = [model_path, tmpf]
                collect_output[cmd, gene_name, foutput]
                
    print(f"Results saved to {foutput}")
    
if __name__ == "__main__":
    main()