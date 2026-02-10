#!/usr/bin/env python3

import sys
import argparse

def parse_region(region_str):
    """Parse chrom:start-end into (chrom, start, end)."""
    chrom, span = region_str.split(":")
    beg, end = map(int, span.split("-"))
    return chrom, beg, end

def load_fasta_regions(fasta_file):
    """Extract regions from FASTA headers."""
    name = ""
    regions_by_chrom = {}
    with open(fasta_file) as f:
        for line in f:
            if not line.startswith(">"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue
            name = parts[0]
            region_str = parts[1]  # chrom:start-end
            chrom, beg, end = parse_region(region_str)
            regions_by_chrom.setdefault(chrom, []).append((beg, end))

    return name[1:].split("_")[0],regions_by_chrom

def overlaps(a_start, a_end, b_start, b_end):
    """Check if two intervals overlap (0-based, right-open)."""
    return a_start < b_end and b_start < a_end

def main():
    parser = argparse.ArgumentParser(description="Find overlapping headers between a FASTA and a header file.")
    parser.add_argument("-i", "--input", required=True, help="Header file to filter")
    parser.add_argument("-f", "--fasta", required=True, help="FASTA file with regions to compare")
    parser.add_argument("-o", "--output", required=True, help="Output file with overlapping headers")
    args = parser.parse_args()

    groupname, regions_by_chrom = load_fasta_regions(args.fasta)

    with open(args.input) as hin, open(args.output, "w") as hout:
        for line in hin:
            if not line.strip() or not line.startswith(">"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue
            region_str = parts[1]  # chrom:start-end
            chrom, beg, end = parse_region(region_str)
            if chrom not in regions_by_chrom:
                continue
            for f_beg, f_end in regions_by_chrom[chrom]:
                if overlaps(beg, end, f_beg, f_end):
                    rest = "_".join(line.split("_")[1:])
                    line = f">{groupname}_{rest}"
                    hout.write(line)
                    break

if __name__ == "__main__":
    main()

