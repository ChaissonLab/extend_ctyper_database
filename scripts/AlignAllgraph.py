#!/usr/bin/env python3

from __future__ import annotations

import argparse
import collections as cl
import pathlib
import subprocess
import os
from typing import Dict, List, Tuple
import multiprocessing as mp

from LocalGraphAlign import GetRegions

BLAST_IDENTITY = 90
BLAST_THREADS  = 1
OVERLAP_GAP    = 30000
EXTENDSIZE = 30000

# globals for workers
QUERYSEQ = None
OUTPUT_PATH = None
OUT_LOCK = None


def reverse_complement(seq: str) -> str:
    trans_table = str.maketrans("ATCGatcg", "TAGCtagc")
    return seq.translate(trans_table)[::-1]


def sh(cmd: str, *, check: bool = False) -> None:
    subprocess.run(cmd, shell=True, check=check)
    
    
def read_fasta(path: pathlib.Path) -> Dict[str, str]:
    seqs: Dict[str, str] = {}
    with path.open() as fh:
        head, buf = "", []
        for ln in fh:
            ln = ln.rstrip()
            if ln.startswith(">"):
                if head:
                    seqs[head] = "".join(buf)
                head, buf = ln[1:].split()[0], []
            else:
                buf.append(ln)
        if head:
            seqs[head] = "".join(buf)
    return seqs


def output_regions(results: List[List], output: str):
    # use global lock
    with OUT_LOCK:
        with open(output, mode="a") as w:
            for result in results:
                chrom, strd, start, end = result
                w.write(f"{chrom}\t{strd}\t{start}\t{end}\n")
                
                
def process_block(
    gindex: int,
    hotspot: Dict[str, List[List[int]]],
    sample: str,
    target: str,
    ):
    # use globals in worker
    query = QUERYSEQ
    output = OUTPUT_PATH
    
    folder = os.path.dirname(target)
    os.makedirs(os.path.join(folder, "samples"), exist_ok=True)
    hotspot_fa = f"{folder}/samples/{sample}_hotspot.txt.fa"
    
    # your hotspots is already grouped by contig: {contig: [[beg,end], ...]}
    regions = hotspot  # already extended in main
    
    allresults = []
    index = 1
    
    for contig, regionlist in regions.items():
        seq_ref = query[contig]
        for beg, end in regionlist:
            seq = seq_ref[beg:end]
            
            # write query fasta for this hotspot
            with open(hotspot_fa, "w") as out:
                out.write(f">{sample}_{index}\t{contig}:{beg}-{end}+\n")
                out.write(seq + "\n")
            index += 1
            
            lowercase = sum(1 for c in seq if c.islower())
            
            # NOTE: adjust args to GetRegions to match your actual signature
            # I assume: GetRegions(input_fa, graph_target, out_tmp, lowercase)
            # we can write to a temp file and read back, but since your original
            # code returned results, let's keep that pattern:
            results = GetRegions(hotspot_fa, target, None, lowercase)
            
            allresults.extend([[contig] + x for x in results])
            
    if allresults:
        output_regions(allresults, output)
        
        
def cli() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Hotspot → local graph → regions")
    p.add_argument("-i", "--input",   required=True, help="hotspot txt (per locus)")
    p.add_argument("-q", "--query",   required=True, help="reference FASTA")
    p.add_argument("-t", "--target",  required=True, help="file listing graph targets, one per line")
    p.add_argument("-s", "--sample",  required=True, help="sample name")
    p.add_argument("-o", "--output",  required=True, help="output BED")
    p.add_argument("-p", "--procs",   type=int, default=4, help="processes")
    return p.parse_args()


def main() -> None:
    global QUERYSEQ, OUTPUT_PATH, OUT_LOCK
    
    args = cli()
    
    # 1) targets
    targets = []
    with open(args.target) as f:
        for line in f:
            line = line.strip()
            if line:
                targets.append(line)
                
    # 2) reference
    QUERYSEQ = read_fasta(pathlib.Path(args.query))
    
    # 3) hotspots input -> build structure:
    #    hotspots[gindex][contig] = [[beg, end], ...]
    hotspots = cl.defaultdict(lambda: cl.defaultdict(list))
    with open(args.input) as f:
        for line in f:
            parts = line.strip().split()
            # expected: chrom ... beg end
            # adjust field indices to your actual format
            contig = parts[0]
            beg0 = int(parts[-2])
            end0 = int(parts[-1])
            
            beg = max(0, beg0 - EXTENDSIZE)
            end = min(len(QUERYSEQ[contig]), end0 + EXTENDSIZE)
            
            gindex = int(parts[2])  # you did this in your code
            hotspots[gindex][contig].append([beg, end])
            
    # 4) prepare output
    OUTPUT_PATH = args.output
    open(OUTPUT_PATH, "w").close()
    
    # 5) shared lock
    OUT_LOCK = mp.Lock()
    
    # 6) run in parallel
    tasks = []
    for gindex, hotspot in hotspots.items():
        # gindex are 1-based, target file is 0-based in your code
        target = targets[gindex - 1]
        tasks.append((gindex, hotspot, args.sample, target))
        
    with mp.Pool(processes=args.procs, initializer=_init_worker,
                initargs=(QUERYSEQ, OUTPUT_PATH, OUT_LOCK)) as pool:
        pool.starmap(process_block, tasks)
    
    
def _init_worker(queryseq, output_path, lock):
    # set globals inside worker
    global QUERYSEQ, OUTPUT_PATH, OUT_LOCK
    QUERYSEQ = queryseq
    OUTPUT_PATH = output_path
    OUT_LOCK = lock
    
    
if __name__ == "__main__":
    main()

