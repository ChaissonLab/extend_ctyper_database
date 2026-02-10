#!/usr/bin/env python3
import argparse
import re

def parse_source_coords(sourcefile):
    """
    Parse source FASTA headers.
    Expected format:
      >seqname chrom:start-end[strand]
    Returns:
      {seqname: (chrom, start, end, strand)}
    """
    coord_map = {}
    with open(sourcefile) as f:
        for line in f:
            if not line.startswith(">"):
                continue
            header = line[1:].strip()
            parts = header.split()
            if len(parts) < 2:
                continue
            seqname, coord_part = parts[0], parts[1]
            match = re.match(r"([^:]+):(\d+)-(\d+)([+-]?)", coord_part)
            if not match:
                raise ValueError(f"[ERROR] Invalid coordinate in source header: {header}")
            chrom, start, end, strand = match.groups()
            strand = strand if strand else "+"
            coord_map[seqname] = (chrom, int(start), int(end), strand)
    return coord_map


def compute_subcoord(src_start, src_end, sub_start, sub_end, strand):
    """
    Compute absolute subregion coordinates given both 0-based inclusive coordinates.
    """
    if strand == "+":
        new_start = src_start + sub_start
        new_end   = src_start + sub_end
    else:
        new_start = src_end - sub_end
        new_end   = src_end - sub_start
    return new_start, new_end, strand


def update_headers(inputfile, source_coords, outputfile):
    with open(inputfile) as fin, open(outputfile, "w") as fout:
        for line in fin:
            if not line.startswith(">"):
                fout.write(line)
                continue

            header = line[1:].strip().rstrip(",")
            # Expect name_start_end
            try:
                name, start, end = header.rsplit("_", 2)
            except ValueError:
                raise ValueError(f"[ERROR] Invalid header format: {header}")

            start, end = int(start), int(end)

            if name not in source_coords:
                raise KeyError(f"[ERROR] {name} not found in source file.")

            chrom, src_start, src_end, strand = source_coords[name]
            new_start, new_end, strand = compute_subcoord(
                src_start, src_end, start, end, strand
            )

            fout.write(f">{name} {chrom}:{new_start}-{new_end}{strand}\n")
    print(f"[INFO] Updated FASTA written to {outputfile}")


def main():
    parser = argparse.ArgumentParser(
        description="Edit FASTA headers based on source coordinates (0-based inclusive)."
    )
    parser.add_argument("-i", "--inputfile", required=True, help="Input FASTA file.")
    parser.add_argument("-s", "--sourcefile", required=True, help="Source FASTA file.")
    parser.add_argument("-o", "--outputfile", required=True, help="Output FASTA file.")
    args = parser.parse_args()

    source_coords = parse_source_coords(args.sourcefile)
    update_headers(args.inputfile, source_coords, args.outputfile)


if __name__ == "__main__":
    main()

