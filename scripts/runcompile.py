#!/usr/bin/env python3

import os
import subprocess
import sys

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
    rerun = 0
    inputfile = sys.argv[1]
    kmerfile = inputfile + "_kmer.list"
    normfile = inputfile + "_norm.gz"
    treefile = inputfile + "_tree.fa"
    newkmerfile = inputfile + "_tree.fa_kmer.list"
    matrixfile = inputfile + "_kmatrix.txt"
    
    print("running: " + inputfile)
    
    folder = os.path.dirname(os.path.abspath(inputfile))
    script_folder = os.path.dirname(os.path.abspath(__file__))
        
    # Step 1: skip if matrix file exists and is large enough
    if not rerun and (os.path.isfile(matrixfile) and os.path.getsize(matrixfile) > 100):
        return
    
    # Step 2: run kmernorm
    if rerun or not (os.path.isfile(normfile) and os.path.getsize(normfile) > 100):
        run_or_fail(f"{script_folder}/kmernorm -i {inputfile} -k {kmerfile} -o {normfile} -w 1 -h 1")
        
    # Step 4: run kmertree
    if rerun or not (os.path.isfile(treefile) and os.path.getsize(treefile) > 100):
        run_or_fail(f"{script_folder}/kmertree -i {inputfile} -n {normfile} -o {treefile}")
        
    # Step 5: run kmerannotate
    if rerun or not (os.path.isfile(newkmerfile) and os.path.getsize(newkmerfile) > 100):
        gaf_file = inputfile + "_lineargraph.gaf"
        if os.path.isfile(gaf_file) and os.path.getsize(gaf_file) > 100:
            run_or_fail(f"python {script_folder}/kmerannotate.py -k {kmerfile} -s {inputfile} -g {inputfile}_graph.FA -a {inputfile}_allgraphalign.out -o {newkmerfile}")
        else:
            run_or_fail(f"python {script_folder}/kmerannotate.py -k {kmerfile} -o {newkmerfile}")
            
    # Step 6: compile matrix
    if rerun or not (os.path.isfile(matrixfile) and os.path.getsize(matrixfile) > 100):
        run_or_fail(f"python {script_folder}/matrixcompile.py -s {treefile} -k {newkmerfile} -o {matrixfile}")
        
if __name__ == "__main__":
    main()
