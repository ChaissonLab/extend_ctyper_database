#!/usr/bin/env python3

import os
import threading
import argparse
import concurrent.futures
import datetime
import collections as cl
import re
import glob
from gfixbreaks import graphDB, assemblysmall, readgraph, qposi_togposi, MergeChunks
global mergesmall
import time


script_dir = os.path.dirname(os.path.abspath(__file__))
cutoffdistance = 30000
REPEAT_CUT=5000
BLAST_IDENTITY = 0.9

def Overlap2(query, ref):
    
    nquery = 2*len(query)
   
    all_coordinates = [x for y in query + ref for x in y[:2]]
    
    sort_index = sorted(range(len(all_coordinates)), key = lambda x: all_coordinates[x])
    
    allgroups = []
    
    number_query_seq = 0
    number_ref_seq = 0
    for index in sort_index:
        
        coordinate = all_coordinates[index]
        
        if index %2 == 0 :
            if index < nquery:
                number_query_seq += 1
                if number_query_seq == 1 and number_ref_seq>0:
                    allgroups.append([coordinate,coordinate])

            else:
                number_ref_seq += 1
                if number_ref_seq == 1 and number_query_seq > 0:
                    allgroups.append([coordinate,coordinate])
                
                
        else:
            if index < nquery:
                number_query_seq -= 1
                if number_query_seq == 0 and number_ref_seq>0:
                    allgroups[-1][1] = coordinate


            else:
                number_ref_seq -= 1
                if number_ref_seq ==0 and number_query_seq > 0:
                    allgroups[-1][1] = coordinate
                
                
                
    return allgroups

def overlapvalid(regions, graphinfo, validregions):
    
    index = 1
    newregions = []
    for locus, regions in regions.items():
        
        for region in regions:
            
            localstart,localend = region[3], region[4]
            
            graphcoordi = qposi_togposi(graphinfo[locus], [localstart,localend-1])
            
            graphspans = cl.defaultdict(list)
            for index in range(len(graphcoordi)//2):
                
                path, start = graphcoordi[2*index][:2]
                path, end = graphcoordi[2*index+1][:2]
                
                graphspans[path].append(sorted([start, end+1]))
                
            totalvalidoverlap = 0
            for path, pathregions in graphspans.items():
                
                validoverlap = Overlap2(pathregions,[validregions[path]])
                totalvalidoverlap += sum([x[1]-x[0] for x in validoverlap])
                
            if totalvalidoverlap > min(1000,0.5*(localend-localstart)):
                newregions.append([f"seq_{index}",locus,localstart,localend,"",""])
                index += 1
                
    return newregions


def clean_nonvalid(breaksoncontigs, validregions):
    
    contigsizes = overlapvalid(breaksoncontigs, validregions)
    
    breaksoncontigs = [y for k,v in breaksoncontigs.items() for y in [[k,x[3],x[4]] for x in v] if len(v)]
    
    iffix = 1
    while iffix:
        breaksoncontigs, iffix = assemblysmall(breaksoncontigs)
    
    return breaksoncontigs


def uniformbreaks(alignout, cachefile, mergeall):
        
    global mergesmall
    mergesmall = 20000

    thegraph = graphDB()
    thegraph = thegraph.load_json(cachefile)
    
    graphinfo, contigspan = readgraph(alignout)
    results = dict()

    for loci in contigspan.keys():
        results[loci] = thegraph.overlap_genes(contigspan[loci], graphinfo[loci])
    
    iffix = 1
    while iffix:
        results, iffix = assemblysmall(results, mergesmall)

    regions = overlapvalid(results,  graphinfo, thegraph.validregions)
        

    #maxsize, regions = MergeChunks(regions, mergeall)

    return regions




def alignment_graph(fasta, graph, alignout, lowercase, nthreads) -> None:
    
    ifrepeat = 1 if lowercase > REPEAT_CUT else 0
    repeatopts = "-dust yes -lcase_masking" if ifrepeat else ""
    
    cmd = f"{script_dir}/KmerStrd -i {fasta} -r {graph} -o {fasta}_ && mv {fasta}_ {fasta}"
    os.system(cmd)
    graph_db = graph + "_db"
    if not os.path.isfile(graph_db + ".ndb"):
        os.system(f"bash {script_dir}/runmakeblastdb -in {graph} -out {graph_db}")

    """Run BLASTN against genome DB and exon DB, appending to *.out."""
    cmd = f"bash {script_dir}/runblastn -task megablast -query {fasta} -db {graph_db} -gapopen 10 -gapextend 2 -word_size 30   -perc_identity {BLAST_IDENTITY} {repeatopts} -evalue 1e-200  -outfmt 17 -out {alignout}  -num_threads {nthreads} -max_target_seqs 100 "
    os.system(cmd)
    
    cmd = f"{script_dir}/KmerStrd -in {fasta} -r {graph} -o {fasta}_"

    if ifrepeat:
        cmd = f"{script_dir}/runwinnowmaplight.sh {fasta} {graph} {nthreads} {BLAST_IDENTITY} 300 >> {alignout}"
        os.system(cmd)

    cmd = f"python {script_dir}/graphcigarlight.py -i {alignout} -q {fasta} -r {graph} -o {alignout}.txt"
    os.system(cmd)




def GraphAlign(fasta, graph, output, lowercase, mergeall):
           
    alignout = output  

    cachefile = graph.replace(".FA","cache.json")

    regions = uniformbreaks(alignout, cachefile,mergeall)
    
    with open(fasta, mode = 'r') as f:
        strd = f.readline().split()[1][-1]

    regions = [[strd,x[2],x[3]] for x in regions]

    return regions

def GetRegions(input,graph,output,lowercase, threads = 1, mergeall = 80000):
    
    alignout = output 
    #t0 = time.time()
    alignment_graph(input, graph, alignout, lowercase, threads)
    #print("alignment_graph:", round(time.time() - t0, 3), "s")

    #t0 = time.time()
    regions = GraphAlign(input, graph, output+".txt", lowercase, mergeall)
    #print("GraphAlign:", round(time.time() - t0, 3), "s")

    return regions

def main(args):

    regions = GetRegions(args.input,args.graph,args.output, args.lowercase, args.threads, args.mergeall)
    print(regions)
    
def run():
    """
        Parse arguments and run
    """
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-i", "--input", help="path to input data file",dest="input", type=str, required=True)
    parser.add_argument("-g", "--graph", help="path to input data file",dest="graph", type=str, required=True)
    parser.add_argument("-o", "--output", help="path to output file", dest="output",type=str, required=True)
    parser.add_argument("-t", "--threads", help="path to output file", dest="threads",type=int, default = 1)
    parser.add_argument("-l", "--lowcase", help="path to output file", dest="lowercase",type=int, default = 0)   
    parser.add_argument("-m", "--mergeall", help="path to output file", dest="mergeall",type=int, default = 80000)
    
    parser.set_defaults(func=main)
    args = parser.parse_args()
    args.func(args)
    
    
if __name__ == "__main__":
    run()
