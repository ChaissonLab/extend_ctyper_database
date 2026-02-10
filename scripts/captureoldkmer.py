#!/usr/bin/env python3


import os
import subprocess
import sys

script_folder = os.path.dirname(os.path.abspath(__file__))

def main():
    
    file1 = sys.argv[1]
    file2pref = sys.argv[2]
    file2 = "../seventhrunfix2/groups/{}/{}.fasta".format(file2pref, file2pref)
    file2kmer = file2 + "_kmer.list"
    folder = "../seventhrunfix2/groups/{}/".format(file2pref)
    
    file3= file1.split("_partitions")[0]
    file3kmer = file3 + "_kmer.list"
    
    os.system("kmer_convertor -i {} -t {} -o {}/old.counts".format(file3, file3kmer, folder))
    os.system("kmer_convertor -i {} -t {} -o {}/new.counts".format(file2, file3kmer, folder))
    
    old_counts_file = os.path.join(folder, "old.counts")
    new_counts_file = os.path.join(folder, "new.counts")
    same_counts_file = os.path.join(folder, "same.counts")
    
    # Step 1: find common lines
    with open(old_counts_file) as f:
        f.readline()
        old_lines = set(line.strip() for line in f)
        
    with open(new_counts_file) as f:
        f.readline()
        new_lines = set(line.strip() for line in f)
            
    # Step 2: load file2kmer first-column values
    with open(file2kmer) as f:
        kmer_keys = set(line.strip().split()[0] for line in f if line.strip() and not line.startswith(">"))
    
    same_lines = old_lines & new_lines
    # Step 3: filter same.counts where column 1 not in file2kmer
    with open(same_counts_file, 'w') as out:
        for line in same_lines:
            if line.split()[0] in kmer_keys:
                continue
            out.write(">\n"+line.split()[0]+"\n")
    
    
        
if __name__ == "__main__":
    main()
