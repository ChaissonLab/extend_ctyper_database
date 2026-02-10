#!/usr/bin/env python3

import os
import sys
import argparse
import collections as cl

class error_region:
    def __init__(self, loadfile, extend = 25000):
        
        self.extend = extend
        self.error_regions = cl.defaultdict(list)
        
        if len(loadfile) == 0:
            return
        
        with open(loadfile,mode = 'r') as f:
            
            for line in f:
                
                if line.startswith(">"):
                    
                    region = line.split()[1]
                    contig, locus = region.split(":")
                    if locus.endswith("+") or locus.endswith("-"):
                        locus = locus[:-1]
                    start,end = "-".join(locus.split("-")[:-1]),locus.split("-")[-1]
                    
                    start,end = int(start),int(end)
                    
                    self.error_regions[contig].append([start-extend,end+extend])
    
    def add(self, contig, region):
        
        self.error_regions[contig].append([region[0] - self.extend, region[1] + self.extend])
    
    def findoverlap(self, contig, region):
        
        error_oncontig = self.error_regions[contig]
        
        if len(error_oncontig) == 0:
            return False
        
        region_size = region[1]-region[0]
        for error_region in error_oncontig:
            
            allcordi = error_region + region
            
            overlap =   max(allcordi) - min(allcordi)  -   (error_region[1]-error_region[0]) - region_size
            if overlap <= 0:
                
                return True
            
        return False

def extendfilter(allcontigs,error_regions):
    
    filtered = set()
    for contig, regions in allcontigs.items():
        
        keepfilter = 1
        while keepfilter:
            keepfilter = 0
            for region in regions:
                if error_regions.findoverlap(contig, region[:2]):
                    if region[2] in filtered:
                        continue

                    keepfilter = 1
                    filtered.add(region[2])
                    error_regions.add(contig, region[:2])
                    break
    
    return filtered
    
def main(args):
    
    error_regions = error_region(args.exclude)
    
    
    iffilter_any = 0
    iffilter = 0
    
    allcontigs = cl.defaultdict(list)
    filtered = set()
    with open(args.input,mode = 'r') as r:
        
        for line in r:
            
            if line.startswith(">"):
                
                name,region = line.split()[0], line.split()[1]
                contig, locus = region.split(":")
                if locus.endswith("-") or locus.endswith("+"):
                    locus = locus[:-1]
                start,end = "-".join(locus.split("-")[:-1]),locus.split("-")[-1]
                start,end = int(start),int(end)
                
                allcontigs[contig].append([start,end,name])
                
    filtered = extendfilter(allcontigs,error_regions)
    
    
    
    filterfile = args.input+"_exfilter"
    tempfile = args.input+"_fixtemp"
    if len(filtered):
        
        with open(args.input,mode = 'r') as r, open(tempfile,mode = 'w') as w, open(filterfile,mode = 'w') as w2:
            
            for line in r:
                
                if line.startswith(">"):
                    
                    name= line.split()[0]
                    iffilter = 0
                    if name in filtered:
                        iffilter = 1
                        
                if not iffilter:
                    w.write(line)
                else:
                    w2.write(line)
                    
        os.system("mv {} {}".format(tempfile, args.input))
        
        
def run():
    """
        Parse arguments and run
    """
    parser = argparse.ArgumentParser(description="program determine psuedogene")
    parser.add_argument("-i", "--input", help="path to input data file",dest="input", type=str, required=True)
    parser.add_argument("-e", "--exclude", help="path to input data file",dest="exclude", type=str,default = "")
    parser.set_defaults(func=main)
    args = parser.parse_args()
    args.func(args)
    
    
if __name__ == "__main__":
    run()
    
