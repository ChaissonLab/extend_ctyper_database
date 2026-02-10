#!/usr/bin/env python3
import os
import sys
import argparse
from collections import defaultdict
from multiprocessing import Pool

def process_file(filepath):
    """Extract haplotypes and lines from a single file."""
    results = defaultdict(list)
    print(filepath)
    with open(filepath) as f:
        for line in f:
            if not line.startswith(">"):
                continue
            first_col = line.strip().split("\t")[0]
            parts = first_col.split("_")
            if len(parts) >= 3:
                haplotype = "_".join(parts[1:3])
                results[haplotype].append(line.strip())
    return results

def main(input_list, output_dir, processes):
    os.makedirs(output_dir, exist_ok=True)

    # Read list of files
    with open(input_list) as f:
        files = [line.strip().split("\t")[0] for line in f if line.strip()]
        files = [x for x in files if os.path.isfile(x)]

    combined = defaultdict(list)

    print(len(files))
    # Process files in parallel with multiprocessing
    with Pool(processes=processes) as pool:
        for file_results in pool.imap_unordered(process_file, files):
            for hap, lines in file_results.items():
                combined[hap].extend(lines)

    # Write results
    for hap, lines in combined.items():
        outfile = os.path.join(output_dir, f"{hap}.txt")
        with open(outfile, "w") as out:
            out.write("\n".join(lines) + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Summarize haplotype lines into files by haplotype.")
    parser.add_argument("-i", "--input", required=True, help="Path to all_bfixpartitions.list")
    parser.add_argument("-o", "--output", required=True, help="Output folder for haplotype files")
    parser.add_argument("-t", "--threads", type=int, default=8, help="Number of processes (default: 8)")

    args = parser.parse_args()
    main(args.input, args.output, args.threads)

