#!/usr/bin/env python3

import os
import subprocess
import sys
from pathlib import Path

def run_or_fail(command):
    try:
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Command failed: {command}")
        print(f"[ERROR] Exit code: {e.returncode}")
        if e.returncode == -9:
            print("[OOM] Likely killed by out-of-memory (SIGKILL).")
        sys.exit(1)
        
def main():
    inputfile = sys.argv[1]
    kmerfile = inputfile + "_kmer.list"
    normfile = inputfile + "_norm.gz"
    treefile = inputfile + "_tree.fa"
    newkmerfile = inputfile + "_tree.fa_kmer.list"
    matrixfile = inputfile + "_kmatrix.txt"
    


    print("running: " + inputfile)
    
    folder = os.path.dirname(os.path.abspath(inputfile))
    script_folder = os.path.dirname(os.path.abspath(__file__))
   
    folders = str(folder).split("/")
    folder2 = "/".join(str(folder).split("/")[:-4])
    


    graphfile = "{}/{}{}{}_samples.fasta_fixed.fa_loci.txt.fasta".format("/".join(folders[:-2]),folders[-5].split("_")[0], folders[-3], "_".join(folders[-5].split("_")[1:]))
    graphinfo = "{}/{}{}{}_samples.fasta_fixed.fa_loci.txt".format("/".join(folders[:-2]),folders[-5].split("_")[0], folders[-3], "_".join(folders[-5].split("_")[1:]))
    graphkmer = "{}/{}{}{}_samples.fasta_fixed.fa_qc.fa_kmers.txt".format("/".join(folders[:-2]),folders[-5].split("_")[0], folders[-3], "_".join(folders[-5].split("_")[1:]))
    graphfasta = "{}/{}{}{}_samples.fasta_fixed.fa_qc.fa".format("/".join(folders[:-2]),folders[-5].split("_")[0], folders[-3], "_".join(folders[-5].split("_")[1:])) 
    output = "{}/partitions/".format("/".join(folders[:-2]))
    queryfile = "{}/{}_origin.fa".format(folder2,folders[-5])

    print(f"python {script_folder}/partition.py -k {graphkmer} -i {graphfasta}  -c 1000 -o {output} -q {queryfile}  ")
    run_or_fail(f"python {script_folder}/partition.py -k {graphkmer} -i {graphfasta}  -c 1000 -o {output} -q {queryfile}  ")

            
        
if __name__ == "__main__":
    main()
