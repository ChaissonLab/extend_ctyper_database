#!/usr/bin/env python3

import collections as cl
import os 
import argparse

def reassemble(inputfile, outputfile, normfile):
    
    annotations = dict()
    regions = cl.defaultdict(list)
    contigstrd = cl.defaultdict(int)
    with open(inputfile, mode = 'r') as f:
        for line in f:
            if line.startswith(">"):
                line = line.strip().split()
                name = line[0][1:]
                region = line[1]
                contig, coordi = region.split(":")
                coordi = [int(x) for x in coordi[:-1].split("-")]
                strd = line[1][-1]
                size = coordi[1] - coordi[0]
                regions[contig].append(coordi+[[name]])
                annotations[name] = [size,strd] + line[2:]
                
    regions = {contig: sorted(region, key= lambda x: (x[0],x[1]) if annotations[x[2][0]][1] == "+" else (x[1],x[0]))  for contig,region in regions.items()}
    
    
    ifmerge = 0
    

    new_regions = []
    for contig, regions_oncontig in regions.items():
       
        new_regions_oncontig = []
        lastend = -10000000
        lastsize = 10000000
        for index, region in enumerate(regions_oncontig):
          
            size = region[1] - region[0]
            dis = region[0] - lastend
            if min(size,lastsize) < 10000 and dis < 10000:
                new_regions_oncontig[-1][1] = region[1]
                new_regions_oncontig[-1][2].extend(region[2])
                lastsize += (dis+size)
                ifmerge = 1
            else:
                new_regions_oncontig.append(region)
                lastsize = size 
            lastend = region[1]
        
       
        for i, new_region in enumerate(new_regions_oncontig):
            largest_name = sorted(new_region[-1], key = lambda name: annotations[name], reverse = 1)[0]
            strd = annotations[largest_name][1]
            ifexon = "Exon" if any(x for x in new_region[-1] if annotations[x][2] == "Exon") else "Intron"
            chm13 = ";".join(list(set([ annotations[x][3] for x in new_region[-1]  ])))
            
            new_regions.append([largest_name,contig, strd] + new_region[:2] + [ifexon, chm13])

        
        
    return ifmerge,new_regions


def process_locus2(region, haplopath,outputfile):
    
    newname,  contig, strd,start, end, ifexon, mapinfo = region
   
    if "#" not in contig: 
        haplo = "CHM13_h1" if "NC_0609" in contig else "HG38_h1"
    else:
        haplo = "_h".join(contig.split("#")[:2])
        
    path = haplopath[haplo]
        
    header = ">{}\t{}:{}-{}{}\t{}\t{}".format(newname, contig, start, end, strd, ifexon, mapinfo)
    
    cmd1 = f'echo "{header}" >> {outputfile} && samtools faidx {path} {contig}:{start}-{end} | tail -n +2 >> {outputfile}'
    
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
    
    if 1==1:
        for region in regions:
            # Submit jobs to executor
            #process_locus2( region, haplopath, folder, outputfile)
            process_locus2(region, haplopath,  outputfile)
            # Wait for all threads to complete
            #concurrent.futures.wait(futures)

def main(args):
    
    ifmerge, newregions = reassemble(args.input, args.output, args.query)
   
    if ifmerge:
        
        makenewfasta2(newregions, args.query, args.output)
    
    else:
        
        os.system("ln -f {} {}".format(args.input, args.output))

    
def run():
    """
        Parse arguments and run
    """
    parser = argparse.ArgumentParser(description="program distract overlap genes and alignments on contigs")
    
    parser.add_argument("-i", "--input", help="path to output file", dest="input", type=str,required=True)
    parser.add_argument("-o", "--output", help="path to output file", dest="output", type=str,required=True)
    parser.add_argument("-q", "--query", help="path to output file", dest="query", type=str,required=True)
    
    parser.set_defaults(func=main)
    args = parser.parse_args()
    args.func(args)
    
    
if __name__ == "__main__":
    run()
