#!/usr/bin/env python3
"""
Integrated hotspot pipeline.

Steps per block directory
-------------------------
1. Build hotspot FASTA from   <block>/hotspots/<sample>_hotspot.txt
2. Run BLASTN (genome + exon DB) → PAF → filtered PAF
3. Generate breakpoint table   via getbreakpoints.py
4. Append extracted regions to <block>/<block>_samples.fasta
    (thread-safe with file locking)

Use --touch-align / --touch-region for Snakemake sentinels.
"""

from __future__ import annotations

import argparse
import collections as cl
import pathlib
import subprocess
import sys
import time
import uuid
from typing import Dict, List, Tuple

import pandas as pd

# --------------------------------------------------------------------------- #
# Tunables
# --------------------------------------------------------------------------- #
BLAST_IDENTITY = 90
BLAST_THREADS  = 4
OVERLAP_GAP    = 20_000   # bp allowed between merged loci
SENTINEL_EMPTY = 5        # file size threshold in bytes

# --------------------------------------------------------------------------- #
# Generic helpers
# --------------------------------------------------------------------------- #
def sh(cmd: str, *, check: bool = False) -> None:
    """Run a shell command with I/O passthrough."""
    subprocess.run(cmd, shell=True, check=check)
    
    
def read_fasta(path: pathlib.Path) -> Dict[str, str]:
    """Load a single-line FASTA into {header: sequence}."""
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


def merge_intervals(spans: List[Tuple[int, int]], gap: int = OVERLAP_GAP) -> List[Tuple[int, int]]:
    """Greedy merge of (start,end) pairs allowing ≤ *gap* bp separation."""
    if not spans:
        return []
    
    spans.sort()
    merged = [list(spans[0])]
    for s, e in spans[1:]:
        last = merged[-1]
        if s - last[1] <= gap:
            last[1] = max(last[1], e)
        else:
            merged.append([s, e])
            
    return [(max(0, s - gap // 2), e + gap // 2) for s, e in merged]


def run_blast(fasta: pathlib.Path, db: pathlib.Path) -> None:
    """Run BLASTN against genome DB and exon DB, appending to *.out."""
    out = f"{fasta}_blast.out"
    sh(
        f"blastn -query {fasta} -db {db}.fa_db "
        f"-perc_identity {BLAST_IDENTITY} -dust yes -lcase_masking "
        f"-outfmt 17 -num_threads {BLAST_THREADS} -out {out}"
    )
    sh(
        f"blastn -query {fasta} -db {db}_exons.fa_db "
        f"-perc_identity {BLAST_IDENTITY} -outfmt 17 -num_threads {BLAST_THREADS} "
        f"-out - >> {out}"
    )
    
# --------------------------------------------------------------------------- #
# Locked-file helper
# --------------------------------------------------------------------------- #
class LockedFile:
    """Simple process-level file lock for concurrent appends."""
    def __init__(self, filepath: pathlib.Path, mode: str = "a", max_age: int = 60):
        self.filepath = filepath
        self.mode     = mode
        self.max_age  = max_age
        self.lockfile = filepath.with_suffix(f".lock_{uuid.uuid4().hex[:8]}")
        self._fh      = None
        
    def _cleanup(self) -> None:
        now = time.time()
        for lf in self.filepath.parent.glob(f"{self.filepath.name}.lock_*"):
            try:
                if now - lf.stat().st_mtime > self.max_age:
                    lf.unlink(missing_ok=True)
            except FileNotFoundError:
                pass
                
    def __enter__(self):
        while True:
            self._cleanup()
            try:
                self.lockfile.write_text(str(time.time()))
                break
            except FileExistsError:
                time.sleep(1)
        self._fh = self.filepath.open(self.mode)
        return self._fh
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        if self._fh:
            self._fh.close()
        self.lockfile.unlink(missing_ok=True)
        
# --------------------------------------------------------------------------- #
# Pipeline helpers
# --------------------------------------------------------------------------- #
def build_hotspot_fasta(bed: pathlib.Path, fasta_out: pathlib.Path, query: Dict[str, str]) -> None:
    """Convert <sample>_hotspot.txt → FASTA of merged loci."""
    try:
        df = pd.read_csv(bed, sep="\t", header=None)
    except FileNotFoundError:
        print(f"[WARN] hotspot file missing: {bed}", file=sys.stderr)
        return
    
    groups: Dict[str, List[int]] = cl.defaultdict(list)
    for idx, contig in enumerate(df.iloc[:, 0].str.split().str[0]):
        groups[contig].append(idx)
        
    with fasta_out.open("w") as out:
        for contig, idxs in groups.items():
            spans = [tuple(map(int, x)) for x in df.iloc[idxs, -2:].values]
            for beg, end in merge_intervals(spans):
                out.write(f">{contig}_{beg}_{end}\t{contig}:{beg}-{end}\n")
                out.write(query[contig][beg:end] + "\n")
                
                
def run_breakpoints_and_append(
    txt_out: pathlib.Path,
    fasta_fn: pathlib.Path,
    sample: str,
    query: Dict[str, str],
    ) -> None:
    """Append breakpoint regions into *fasta_fn* (thread-safe)."""
    if not txt_out.exists() or txt_out.stat().st_size < SENTINEL_EMPTY:
        return
    
    with txt_out.open() as src, LockedFile(fasta_fn, "a") as dst:
        counter = 1
        prefix  = fasta_fn.stem.split("_")[0]
        for ln in src:
            parts = ln.split()
            if len(parts) < 4:
                continue                    # malformed line
            title, chrom, beg, end = parts[:4]
            exon  = parts[6] if len(parts) > 6 else "."
            header = f"{prefix}_{sample}_{counter}\t{chrom}:{beg}-{end}\t{exon}\t{title}"
            seq    = query[chrom][int(beg):int(end)]
            dst.write(f">{header}\n{seq}\n")
            counter += 1
            
            
def process_block(
    block: pathlib.Path,
    sample: str,
    query: Dict[str, str],
    undo: bool,
    script_dir: pathlib.Path,
    ) -> None:
    hotspot_txt = block / "hotspots" / f"{sample}_hotspot.txt"
    hotspot_fa  = hotspot_txt.with_suffix(".txt.fa")
    db_prefix   = block / block.name
    
    paf        = hotspot_fa.with_suffix(".fa.paf")
    paf_filt   = hotspot_fa.with_suffix(".fa.paf_filter.PAF")
    title_txt  = hotspot_fa.with_suffix(".fa_title.txt")
    blast_out  = hotspot_fa.with_suffix(".fa_blast.out")
    txt_out    = paf_filt.with_suffix(".PAF.txt")
    query_fn   = block / f"{block.name}_origin.fa"
    
    # 1  Build FASTA
    build_hotspot_fasta(hotspot_txt, hotspot_fa, query)
    
    # OPTIONAL: skip empty blocks to save time (remove if you prefer)
    # if not hotspot_fa.exists() or hotspot_fa.stat().st_size < SENTINEL_EMPTY:
    #     return
    
    # 2  BLAST → PAF → filter
    if undo or not (blast_out.exists() and blast_out.stat().st_size > SENTINEL_EMPTY):
        run_blast(hotspot_fa, db_prefix)
        
    sh(f"python {script_dir}/blastntopaf.py -i {blast_out} -q {hotspot_fa} > {paf}")
    sh(f"python {script_dir}/filter_alignment.py -i {paf} -r {db_prefix}.fa -o {paf_filt}")
    sh(f"grep '^>' {hotspot_fa} > {title_txt} || true")
    
    # 3  Breakpoints
    if paf_filt.exists() and paf_filt.stat().st_size >= SENTINEL_EMPTY:
        sh(f"python {script_dir}/getbreakpoints.py -r {query_fn} -i {paf_filt} -o {txt_out}")
        
    # 4  Append regions to per-block FASTA
    samples_fa = block / f"{block.name.replace('_', '')}_samples.fasta"
    run_breakpoints_and_append(txt_out, samples_fa, sample, query)
    
    hotspot_fa.unlink(missing_ok=True)   # tidy
    
    
# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def cli() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Hotspot → BLAST → region-FASTA pipeline")
    p.add_argument("-i", "--input",   required=True, help="Root folder with block dirs")
    p.add_argument("-q", "--query",   required=True, help="Reference FASTA for this sample")
    p.add_argument("-n", "--sample",  required=True, help="Sample prefix (matches hotspot files)")
    p.add_argument("-u", "--undo", action="store_true",
                    help="Force re-BLAST even if previous results exist")
    p.add_argument("--touch-align",  help="Sentinel file to touch after align phase")
    p.add_argument("--touch-region", help="Sentinel file to touch after region phase")
    return p.parse_args()


def main() -> None:
    args       = cli()
    root       = pathlib.Path(args.input).expanduser().resolve()
    script_dir = pathlib.Path(__file__).resolve().parent
    query_seqs = read_fasta(pathlib.Path(args.query))
    
    for blk in root.iterdir():
        if blk.is_dir():
            process_block(blk, args.sample, query_seqs, args.undo, script_dir)
            
    if args.touch_align:
        pathlib.Path(args.touch_align).touch()
    if args.touch_region:
        pathlib.Path(args.touch_region).touch()
        
        
if __name__ == "__main__":
    main()
    
