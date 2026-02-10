#!/usr/bin/env python3


import os
import subprocess
import sys

script_folder = os.path.dirname(os.path.abspath(__file__))

def main():
   
    rerun = 1
    inputfile = sys.argv[1]
    kmerfile = inputfile + "_kmers.list"
    normfile = inputfile+"_norms.gz"
    treefile = inputfile+"_treeorders.fa"
    newkmerfile = inputfile+"_treeorders.fa_kmer.list"
    matrixfile = inputfile + "_kmatrix.txt"
    
    print ("running: "+ inputfile) 
    folder = os.path.dirname(os.path.abspath(inputfile))
    extrakmer = folder+"/same.list"
    exclude_file = inputfile + "_kmers.txt_exclude.txt"
    kmer_list_file = inputfile + "_kmer.list"
   
    if not rerun and ( (os.path.isfile(matrixfile) and os.path.getsize(matrixfile) > 100)):
        return 
    os.system(f"grep -F -v -f {exclude_file} {kmer_list_file} | awk '{{ print \">\\n\" $0 }}' > {kmerfile}")
    
    if os.path.isfile(extrakmer):
        os.system(f"grep -F -v -f {exclude_file} {extrakmer} | awk '{{ print \">\\n\" $0 }}' >> {kmerfile}")
    
    if rerun or ( not (os.path.isfile(normfile) and os.path.getsize(normfile) > 100)):
        os.system("{}/kmernorm -i {} -k {} -o {} -w 1 -h 1".format(script_folder, inputfile, kmerfile, normfile))
    if rerun or (not (os.path.isfile(treefile) and os.path.getsize(treefile) > 100)):
        os.system("python {}/kmertree.py -i {} -k {} -n {} -o {}".format(script_folder, inputfile, kmerfile, normfile, treefile))
    if rerun or (not (os.path.isfile(newkmerfile) and os.path.getsize(newkmerfile) > 100)):
        if (os.path.isfile(inputfile + "_lineargraph.gaf") and os.path.getsize(inputfile + "_lineargraph.gaf") > 100):
            os.system("python {}/kmerannotate.py -k {} -s {} -g {}_graph.FA -a {}_allgraphalign.out -o {}".format(script_folder,kmerfile,inputfile,inputfile,inputfile,newkmerfile))
        else:
            os.system("python {}/kmerannotate.py -k {} -o {}".format(script_folder,kmerfile,inputfile,inputfile,inputfile,newkmerfile))
    if rerun or (not (os.path.isfile(matrixfile) and os.path.getsize(matrixfile) > 100)):
        os.system("python {}/matrixcompile.py -s {}  -k {} -o {} ".format(script_folder,treefile,newkmerfile, matrixfile))
        
if __name__ == "__main__":
    main()
