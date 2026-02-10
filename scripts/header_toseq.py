#!/usr/bin/env python3
import argparse
import subprocess
import sys
from pathlib import Path

def load_query_paths(pathfile):
    """Load haplotype → fasta path mappings."""
    qpaths = {}
    with open(pathfile) as f:
        for line in f:
            if not line.strip():
                continue
            hap, fasta = line.strip().split()[:2]
            qpaths[hap] = fasta
    return qpaths


def parse_header(line):
    """
    Parse header line of form:
      >seq1 contig:start-end+
      >seq2 contig:start-end-
      >seq3 contig:start-end
    """
    header = line.strip().lstrip(">")
    parts = header.split()
    if len(parts) < 2:
        raise ValueError(f"Invalid header: {line}")

    coord_part = parts[1]

    # Strand is + (default) or -
    strand = "+"
    if coord_part.endswith("+"):
        coord_part = coord_part[:-1]
        strand = "+"
    elif coord_part.endswith("-"):
        coord_part = coord_part[:-1]
        strand = "-"

    contig, coords = coord_part.split(":")
    start, end = coords.split("-")
    start, end = int(start), int(end)

    # Determine haplotype name
    if "#" in contig:
        haploname = "_h".join(contig.split("#")[:2])
    elif contig.startswith("NC_0609"):
        haploname = "CHM13_h1"
    else:
        haploname = "HG38_h1"

    return header, contig, start, end, strand, haploname


def fetch_sequence(fasta, contig, start, end):
    """Fetch sequence with samtools faidx."""
    region = f"{contig}:{start-1}-{end}"  # samtools expects 1-based, so start-1
    cmd = ["samtools", "faidx", fasta, region]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"samtools faidx failed for {region}")
    seq_lines = result.stdout.splitlines()[1:]  # skip header
    return "".join(seq_lines)


def main():
    parser = argparse.ArgumentParser(description="Convert FASTA headers to sequences using samtools faidx")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA header file")
    parser.add_argument("-q", "--query", required=True, help="query_pathes.txt (haplotype fastapath)")
    parser.add_argument("-o", "--output", required=True, help="Output FASTA file")
    args = parser.parse_args()

    qpaths = load_query_paths(args.query)

    with open(args.output, "w") as out, open(args.input) as f:
        for line in f:
            if not line.startswith(">"):
                continue

            try:
                header, contig, start, end, strand, haploname = parse_header(line)
            except Exception as e:
                print(f"⚠️ Skipping malformed header: {line.strip()} ({e})", file=sys.stderr)
                continue

            fasta_path = qpaths.get(haploname)
            if not fasta_path:
                print(f"⚠️ No FASTA found for {haploname}, skipping", file=sys.stderr)
                continue

            try:
                seq = fetch_sequence(fasta_path, contig, start, end)
                out.write(f">{header}\n{seq}\n")
            except Exception as e:
                print(f"❌ Failed fetching {contig}:{start}-{end} from {haploname} ({e})", file=sys.stderr)
                continue


if __name__ == "__main__":
    main()

