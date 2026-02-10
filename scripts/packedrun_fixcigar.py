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

    graphfile = "{}/{}{}{}_samples.fasta_fixed.fa_loci.txt.fasta".format("/".join(folders[:-2]),folders[-5].split("_")[0], folders[-3], "_".join(folders[-5].split("_")[1:]))
    graphinfo = "{}/{}{}{}_samples.fasta_fixed.fa_loci.txt".format("/".join(folders[:-2]),folders[-5].split("_")[0], folders[-3], "_".join(folders[-5].split("_")[1:]))


    run_or_fail(f"python {script_folder}/getLocalCigar.py -i {treefile} -m {graphfile} -a {graphfile}_allgraphalign.out -o {treefile}_allgraphalign.out")

            
        
if __name__ == "__main__":
    main()
