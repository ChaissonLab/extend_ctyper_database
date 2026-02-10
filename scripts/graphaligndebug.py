#!/usr/bin/env python3

import multiprocessing as mul
import argparse
import os
import subprocess
import datetime
from graphfilesplit import filesplit
from graphinserts import insertdistract
from graphcigar import graphcigar
from graphcleanrepeats import cleanrepeats
from graphtolinear import graphtolinear
import re

script_folder = os.path.dirname(os.path.abspath(__file__))
lock1 = mul.Lock()

ifblast_g = 1


def checkmask(fasta_file):
    
    bash_cmd = f'''
    grep -v '^>' "{fasta_file}" | tr -d '\\n' | tr -cd 'a-z' | wc -c
    '''
    result = subprocess.run(["bash", "-c", bash_cmd], capture_output=True, text=True)
    lowercase_count = int(result.stdout.strip())
    
    return lowercase_count

def fileordered(allfiles, priors=None):
    
    def file_has_Nn(fa_path):
            with open(fa_path, "r") as f:
                for line in f:
                    if line.startswith('>'):   # header
                        continue
                    if ('N' in line) or ('n' in line):
                        return True
            return False
    
    def header_has_chr_colon(fa_path):
        _chr_pat = re.compile(r"^chr(?:[0-9]+|X|Y):")
        with open(fa_path, "r") as f:
            return _chr_pat.match(f.readline().strip()) is not None
        
    priors = priors or set()
    
    def score(path):
        fn = os.path.basename(path)
        stem = fn[:-3] if fn.endswith(".fa") else os.path.splitext(fn)[0]
        sz = os.path.getsize(path)
        
        if stem in priors:
            return 10e20 + sz
        elif ("NC_0609" in fn) or ("CHM13_" in fn):
            return 5e20 + sz
        elif "HG002_" in fn:
            return 10e17 + sz
        if file_has_Nn(path):
            sz = sz - 10e20
        if "HG38alt" in fn:
            return 10e18 + sz
        elif (("_chr" in fn) or ("HG38_" in fn)):
            return (10e19 + sz) if header_has_chr_colon(path) else (10e18 + sz)
        
        else:
            return sz
        
    return sorted(allfiles, key=score, reverse=True)

def selfalign(input,output,nthreads, ifblast = 0):
    
    if ifblast:
        dbcmd = "bash {}/runmakeblastdb -in {} -out {}_db ".format(script_folder,input, input)
        os.system(dbcmd)
        cmd = "bash {}/runblastn  -task megablast  -query {} -db {}_db -gapopen 10 -gapextend 2 -word_size 30  -perc_identity 95 -dust yes -lcase_masking -evalue 1e-200 -outfmt 17 -out {}_selfblast.out  -num_threads {} -max_target_seqs 100 ".format( script_folder,input , input, input,nthreads)
        os.system(cmd)
        totalclean = cleanrepeats(input, input+"_", nthreads)
        
        cmd = "{}/runwinnowmap.sh {} {} {} 0.95 150  > {}_selfblast.out".format(script_folder,  input+"_",input+"_", nthreads, input+"_")
        os.system(cmd)
        totalclean += cleanrepeats(input+"_", output, nthreads)
    else:  
        cmd = "{}/runwinnowmap.sh {} {} {} 0.95 150 > {}_selfblast.out".format(script_folder,  input,input, nthreads, input)
        os.system(cmd)
        
        totalclean = cleanrepeats(input, output, nthreads)
        
    return totalclean

def callblastn(query, graphfile,  output,nthreads, simi = 0.90, ifrepeat = 1):
    
    if ifrepeat:
        repeatopts = "-dust yes -lcase_masking"
    else:
        repeatopts = ""
        
    cmd = "bash {}/runblastn -task megablast -query {} -db {} -gapopen 10 -gapextend 2 -word_size 30   -perc_identity {} {} -evalue 1e-200  -outfmt 17 -out {}  -num_threads {} -max_target_seqs 100 ".format( script_folder, query,graphfile+"_db" , int(simi * 100),repeatopts,output,nthreads)
    os.system(cmd)
    
    
    
def callalign(query, graphfile,  output,nthreads, simi = 0.90,  ifhighqual = 1, ifrepeat = 1):
    if ifhighqual:
        callblastn(query, graphfile,  output,nthreads, simi, ifrepeat)
        
        
        #if ifrepeat:
            #repeatopts = "-dust yes -lcase_masking"
        #else:
            #repeatopts = ""
        
        #cmd = "bash {}/runblastn -task megablast -query {} -db {} -gapopen 10 -gapextend 2 -word_size 30   -perc_identity {} {} -evalue 1e-200  -outfmt 17 -out {}  -num_threads {} -max_target_seqs 100 ".format( script_folder, query,graphfile+"_db" , int(simi * 100),repeatopts,output,nthreads)
        #os.system(cmd)
        
        if (not os.path.exists(output) or os.path.getsize(output) <= 10):
            
            print(f"WARNING: alignment fail :{output}, realigning\n")
            
            callblastn(query, graphfile,  output,nthreads, simi, ifrepeat)
            
            #fallback_cmd = "bash {}/runblastn -task megablast -query {} -db {} -gapopen 10 -gapextend 2 -word_size 30 -perc_identity {} {} -evalue 1e-200 -outfmt 17 -out {} -num_threads {} -max_target_seqs 100".format(script_folder, query, graphfile+"_db", int(simi * 100), repeatopts, output, nthreads)
            
            #os.system(fallback_cmd)
            
            
        if ifrepeat:
            
            cmd = "{}/runwinnowmap.sh {} {} {} {} 150 > {}".format(script_folder,  query,graphfile, nthreads, simi, output+"_wm")
            os.system(cmd)
            if (not os.path.exists(output+"_wm") or os.path.getsize(output+"_wm") <= 10):
                os.system(cmd)
            os.system(f" ( cat  {output}_wm  >> {output} &&  rm {output}_wm  ) || true ")
            
    else:
        #cmd = "{}/runfast.sh {} {} {} 0.95 300 -c > {}".format(script_folder,  query,dbpath, nthreads,output)
        cmd = "{}/runwinnowmap.sh {} {} {} {} 150 > {}".format(script_folder,  query,graphfile, nthreads, simi, output)
        os.system(cmd)

def run(args):
    
    dbcmd = "bash {}/runmakeblastdb -in {} -out {} ".format(script_folder, args.db, args.db+"_db")
    os.system(dbcmd)
        
    callalign(args.input, args.db, args.output, 4)
    
    queryname, fullpath, fullcigar, pathranges,qranges = graphcigar(args.db, args.input , args.output)

    print("\t".join([queryname, fullpath, fullcigar, pathranges,qranges]))

def main():
    """
            Parse arguments and run
    """
    parser = argparse.ArgumentParser(description="program distract overlap genes and alignments on contigs")
    
    parser.add_argument("-i", "--input", help="path to output file", dest="input", type=str,required=True)
    parser.add_argument("-o", "--output", help="if generate fine graph", dest="output", type=str,required=True)
    parser.add_argument("-d", "--db", help="if align to graph", dest="db", type=str,required=True)
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()
