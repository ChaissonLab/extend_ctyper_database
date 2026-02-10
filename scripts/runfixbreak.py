#!/usr/bin/env python3


import argparse
import collections as cl
import shutil
import os
from pathlib import Path
import re
import subprocess

script_folder = Path(__file__).resolve().parent         
            
def main(args):
    inputfile = Path(args.input)
    outputfile = Path(args.output)
    kmerfile = Path(args.kmer)
    normfile = Path(args.norm)
    nthreads = args.nthreads
    wkfolder = inputfile.parent
    query = args.query
        
    outputstr = subprocess.check_output(f"bash {script_folder}/checkmask.sh {inputfile}  ", shell=True, text=True).strip()
    totalchar, muskchar = [int(x) for x in outputstr.split()[:2]]
    if muskchar > 0.5 * totalchar:
        return 
    
    #os.system(f"python {script_folder}/fixbreaks.py -i {inputfile} -n {normfile} ")  
    #fixlog = Path(str(inputfile)+"_fixlog.txt")
    #if fixlog.is_file() and fixlog.stat().st_size > 1:
    if 1==1:
        tempfolder = args.input+"_breakfix/"
        os.makedirs(tempfolder,exist_ok=True)
        os.system(f"python {script_folder}/gfixbreaks.py -i {inputfile} -o {tempfolder}/fix.fa -t {nthreads} -q {query} ") 
        os.system(f"mv {tempfolder}/fix.fa {outputfile}") 
        os.system(f"rm -rf {tempfolder} || true")
    
            
def run():
    parser = argparse.ArgumentParser(description="Program to determine pseudogene")
    parser.add_argument("-i", "--input", required=True, type=str)
    parser.add_argument("-o", "--output", type=str, default="")
    parser.add_argument("-k", "--kmer", type=str, default="")
    parser.add_argument("-n", "--norm", type=str, default="")
    parser.add_argument("-t", "--nthreads", type=int, default=1)
    parser.add_argument("-q", "--query", type=str, default="")
    parser.set_defaults(func=main)
    args = parser.parse_args()
    args.func(args)
    
if __name__ == "__main__":
    run()
