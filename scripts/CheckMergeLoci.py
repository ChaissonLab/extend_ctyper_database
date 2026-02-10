#!/usr/bin/env python3

import os
import threading
import argparse
import concurrent.futures
import datetime
import collections as cl
import re
import numpy as np
import glob


def process_locus2(region, haplopath, outputfile):
    
    newname,  contig, start, end, ifexon, mapinfo = region

    if "#" not in contig: 
        haplo = "CHM13_h1" if "NC_0609" in contig else "HG38_h1"
    else:
        haplo = "_h".join(contig.split("#")[:2])
        
    path = haplopath[haplo]

    start = max(0,int(start))
    header = ">{}\t{}:{}-{}\t{}\t{}".format(newname, contig, start, end, ifexon, mapinfo)

    cmd1 = f'echo "{header}" >> {outputfile} && samtools faidx {path} {contig}:{start+1}-{end} | tail -n +2 >> {outputfile}'
    os.system(cmd1)

    return outputfile


def makenewfasta2(regions, queryfile, outputfile):
    
    haplopath = dict()
    with open(queryfile, mode = 'r') as f:
        for line in f:
            line = line.strip().split()
            haplopath[line[0]] = line[1]
                
    os.system("rm {} || true ".format(outputfile))

    alloutputs = []

    #with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
                    #futures = []
    for region in regions:
            # Submit jobs to executor
            process_locus2( region, haplopath, outputfile)
            #futures.append(executor.submit(process_locus2, region, haplopath, folder, outputfile))
            #alloutputs.append(os.path.join(folder, f"{region[0]}.fa"))
            # Wait for all threads to complete
            #concurrent.futures.wait(futures)
        

def MergeChunks(regions, mergedistance):
    
    regions_sort = sorted(regions, key = lambda x: (x[1], x[2]))
    
    lastcontig = ""
    lastend = -mergedistance - 1
    
    ifmerge = 0
    merged_regions = []
    current_region = []
    for region in regions_sort:
        
        newname,  contig, start, end, ifexon = region[:5]
        
        if contig != lastcontig:
            lastend = -mergedistance - 1
       
        if start - lastend < mergedistance:
            current_region[3] = end
            current_region[4] = "Exon" if ifexon == "Exon" else current_region[4]
            ifmerge = 1
        else:
            merged_regions.append(current_region)
            current_region = region
     
        lastend = end
        lastcontig = contig



    merged_regions = merged_regions[1:]
    merged_regions.append(current_region)
    maxsize = max([x[3]-x[2] for x in merged_regions])
    if maxsize < mergedistance:
        regions = merged_regions
    else:
        ifmerge = 0
    
    return ifmerge, regions

def readfile(inputfile):

    regions = []
    with open(inputfile, mode = 'r') as f:
        
        for line in f:
            
            if line.startswith(">"):
                line = line[1:].split()
                name, locus, ifexon = line[:3]
                    
                strd = "+"
                contig, coordi = locus.split(":")
                if coordi[-1] in ["+","-"]:
                        coordi = coordi[:-1]
                        strd = locus[-1]
                    
                coordi = coordi.split("-")
                start, end = int(coordi[0]), int(coordi[1])
                region = [name, contig, start, end,ifexon] +  line[3:]
                regions.append(region)
    return regions

def main(args):
    
    regions = readfile(args.input)

    ifmerge, regions = MergeChunks(regions, args.mergeall)

    if ifmerge:
        
        makenewfasta2(regions, args.query, args.output)
    

def run():
    """
            Parse arguments and run
    """
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-i", "--input", help="path to input data file",dest="input", type=str, required=True)
    parser.add_argument("-o", "--output", help="path to output file", dest="output",type=str, required=True)
    parser.add_argument("-q", "--query", help="path to output file", dest="query",type=str, default = True)
    parser.add_argument("-m", "--mergeall", help="path to output file", dest="mergeall",type=int, default = 100000)
    
    parser.set_defaults(func=main)
    args = parser.parse_args()
    args.func(args)
    
    
if __name__ == "__main__":
    run()
