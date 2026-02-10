#!/usr/bin/env python3

import argparse
from Bio import SeqIO
import pybedtools

def main():
    parser = argparse.ArgumentParser(
        description="Mask FASTA to lowercase, then unmask (uppercase) regions defined in BED."
    )
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    parser.add_argument("-b", "--bed", required=True, help="BED file with regions to unmask")
    parser.add_argument("-o", "--output", required=True, help="Output FASTA file")
    args = parser.parse_args()

    # Load BED regions
    regions = pybedtools.BedTool(args.bed)

    with open(args.output, "w") as out:
        for record in SeqIO.parse(args.input, "fasta"):
            seq = list(str(record.seq).lower())  # mask whole sequence

            # Apply uppercase for defined regions
            for interval in regions.filter(lambda x: x.chrom == record.id):
                start, end = int(interval.start), int(interval.end)
                for i in range(start, end):
                    if i < len(seq):  # avoid out-of-range
                        seq[i] = seq[i].upper()

            record.seq = type(record.seq)("".join(seq))
            SeqIO.write(record, out, "fasta")

if __name__ == "__main__":
    main()

