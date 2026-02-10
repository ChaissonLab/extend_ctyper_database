#!/usr/bin/env python3

import os
import argparse
import collections as cl
from pathlib import Path

nametochr = {
    'NC_060925.1': 'chr1', 'NC_060926.1': 'chr2', 'NC_060927.1': 'chr3', 'NC_060928.1': 'chr4', 'NC_060929.1': 'chr5',
    'NC_060930.1': 'chr6', 'NC_060931.1': 'chr7', 'NC_060932.1': 'chr8', 'NC_060933.1': 'chr9', 'NC_060934.1': 'chr10',
    'NC_060935.1': 'chr11', 'NC_060936.1': 'chr12', 'NC_060937.1': 'chr13', 'NC_060938.1': 'chr14', 'NC_060939.1': 'chr15',
    'NC_060940.1': 'chr16', 'NC_060941.1': 'chr17', 'NC_060942.1': 'chr18', 'NC_060943.1': 'chr19', 'NC_060944.1': 'chr20',
    'NC_060945.1': 'chr21', 'NC_060946.1': 'chr22', 'NC_060947.1': 'chrX', 'NC_060948.1': 'chrY'
    }


def read_fasta_seq(filepath):
    """Read a single-sequence FASTA into string."""
    seq_lines = []
    with open(filepath) as f:
        for line in f:
            if not line.startswith(">"):
                seq_lines.append(line.strip())
    return "".join(seq_lines)


def write_fasta(filepath, name, seq):
    with open(filepath, "w") as f:
        f.write(f">{name}\n")
        for i in range(0, len(seq), 80):
            f.write(seq[i:i+80] + "\n")


def mask_outside_regions(fasta_path, target_bed, chr_name, extract_start):
    """
    Lowercase all bases except regions from target_bed that overlap the local region.
    `extract_start` = genomic start coordinate of the current extracted region (0-based).
    """
    seq = read_fasta_seq(fasta_path).lower()

    # Parse BED and compute relative coordinates
    with open(target_bed) as f:
        for line in f:
            if line.strip() == "" or line.startswith("#"):
                continue
            chrom, start, end, *rest = line.strip().split()
            if chrom != chr_name:
                continue
            start, end = int(start), int(end)
            rel_start = max(0, start - extract_start)
            rel_end = end - extract_start
            for i in range(rel_start, rel_end):
                if 0 <= i < len(seq):
                    seq = seq[:i] + seq[i].upper() + seq[i+1:]

    # Write masked sequence back
    write_fasta(fasta_path, Path(fasta_path).stem, seq)


def main(args):
    singleref = ""
    queries = {}

    # Read references
    if not args.ref.endswith((".fa", ".fasta")):
        with open(args.ref) as f:
            for line in f:
                line = line.strip().split()
                if len(line) >= 2:
                    queries[line[0]] = line[1]
    else:
        singleref = args.ref

    newfiles = set()
    index = cl.defaultdict(int)

    with open(args.input) as f:
        for line_ in f:
            line = line_.strip().split()
            chr, start, end, name = line[0], int(line[1]), int(line[2]), line[-1]
            index[name] += 1
            outdir = Path(args.output) / name
            os.makedirs(outdir, exist_ok=True)
            fasta_path = outdir / f"{name}.fa"

            if fasta_path.exists():
                fasta_path.unlink()

            title = f">{name}_{index[name]}\t{line_.strip()}"

            # Select query FASTA path
            if not singleref:
                if "#" not in chr:
                    haplo = "CHM13_h1" if chr.startswith("NC_0609") else "HG38_h1"
                else:
                    haplo = "_h".join(chr.split("#")[:2])
                query = queries.get(haplo)
                if not query:
                    raise ValueError(f"No reference found for {haplo}")
            else:
                query = singleref

            # Rename chromosome if needed
            if chr.startswith("NC_0609"):
                chr = nametochr.get(chr, chr)

            # Extract region with samtools
            extract_cmd = f"samtools faidx {query} {chr}:{start+1}-{end} > {fasta_path}"
            os.system(extract_cmd)

            # Apply masking if target is provided
            if args.target:
                mask_outside_regions(fasta_path, args.target, chr, start)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract sequences and optionally mask non-target regions.")
    parser.add_argument("-i", "--input", required=True, help="Input BED-like file with chr start end name")
    parser.add_argument("-r", "--ref", required=True, help="Reference FASTA or mapping file")
    parser.add_argument("-o", "--out", required=True, help="Output directory")
    parser.add_argument("-t", "--target", default="", help="BED file of target regions (mask others)")
    args = parser.parse_args()
    main(args)
