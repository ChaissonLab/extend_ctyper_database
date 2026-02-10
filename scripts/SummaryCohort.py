#!/usr/bin/env python3

import argparse
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(description="Organize cohort data into a summary matrix.")
    parser.add_argument("-l", "--liftover", required=True,
                        help="Large liftover TSV file (e.g., HPRCleave_HG38map.tsv)")
    parser.add_argument("-r", "--bedfile", required=True,
                        help="Reference BED file")
    parser.add_argument("-o", "--output", required=True,
                        help="Output filename prefix")
    # Optional sample list
    parser.add_argument("-s", "--samplelist", required=False,
                        help="Optional file containing list of samples to include (one per line)")
    # Reference label for first column
    parser.add_argument("--ref", required=False, default="HG38",
                        choices=["HG38", "CHM13"],
                        help="Reference tag for first column header (HG38 or CHM13)")
    # Convert list
    parser.add_argument("-c", "--convert", required=False,
                        help=("Optional convert list: 1st col = allele to EXCLUDE as primary, "
                              "2nd col = allele to KEEP as primary. "
                              "Converlist[secondcol].append(firstcol); excludelist.add(firstcol). "
                              "In output, firstcol is never used directly; for any output allele X, "
                              "if X in Converlist, we append all attached alleles with commas."))
    return parser.parse_args()


def get_sample_name(qname: str) -> str:
    """
    Extracts sample name based on logic: '_'.join(Qname.split('_')[1:3])
    """
    parts = qname.split("_")
    return "_".join(parts[1:3])


def main():
    args = parse_args()

    # Decide first-column label based on --ref
    if args.ref.upper() == "CHM13":
        ref_label = "CHM13_h1"
    else:
        ref_label = "HG38_h1"

    # ----------------- load convert list -----------------
    # Converlist[secondcol] = [firstcol1, firstcol2, ...]
    # excludelist = {all firstcol values}
    Converlist = defaultdict(list)
    excludelist = set()

    if args.convert:
        print(f"Reading convert list: {args.convert}...")
        with open(args.convert, "rt") as cf:
            for line in cf:
                line = line.strip()
                if not line:
                    continue
                cols = line.split()
                if len(cols) < 2:
                    continue
                firstcol = cols[0]
                secondcol = cols[1]
                Converlist[secondcol].append(firstcol)
                excludelist.add(firstcol)

        # Deduplicate attached lists (stable)
        for k, vals in Converlist.items():
            Converlist[k] = list(dict.fromkeys(vals))
    # -----------------------------------------------------


    # Rmap: { chrom: { region_name: [start_trimmed, end_trimmed, [qnames]] } }
    Rmap = defaultdict(dict)

    # allsampleindex: { sample_name: column_index }
    allsampleindex = defaultdict(int)

    use_strict_sample_list = False
    current_sample_idx = 0

    # Handle sample list
    if args.samplelist:
        print(f"Reading Sample List: {args.samplelist}...")
        use_strict_sample_list = True
        with open(args.samplelist, "rt") as f:
            for line in f:
                s_name = line.strip()
                if s_name:
                    if s_name not in allsampleindex:
                        allsampleindex[s_name] = current_sample_idx
                        current_sample_idx += 1

    # ----------------- read BED -----------------
    print(f"Reading BED file: {args.bedfile}...")
    with open(args.bedfile, "rt") as f:
        for line in f:
            if not line.strip():
                continue
            row = line.strip().split()
            if len(row) <= 7:
                continue
            if int(row[6]) < 0:
                continue

            chrom = row[0]
            start = int(row[1])
            end   = int(row[2])
            region_name = row[3]

            val1 = int(row[6])
            val2 = int(row[7])

            # Store trimmed coordinates and empty qname list
            Rmap[chrom][region_name] = [start + val1, end - val2, []]

    # ----------------- read Liftover -----------------
    print(f"Reading Liftover file: {args.liftover}...")
    with open(args.liftover, "rt") as f:
        for line in f:
            stripped_line = line.strip()
            if not stripped_line.endswith("Pri"):
                continue

            cols = stripped_line.split("\t")
            if len(cols) < 4:
                continue

            Qname = cols[0]
            Rname = cols[1]
            Rloc  = cols[3]

            Qsample = get_sample_name(Qname)

            # Filtering / sample indexing
            if use_strict_sample_list:
                if Qsample not in allsampleindex:
                    continue
            else:
                if Qsample not in allsampleindex:
                    allsampleindex[Qsample] = current_sample_idx
                    current_sample_idx += 1

            chrom_parts = Rloc.split(":")
            if not chrom_parts:
                continue
            chrom = chrom_parts[0]

            if chrom in Rmap and Rname in Rmap[chrom]:
                Rmap[chrom][Rname][2].append(Qname)

    # ----------------- write output -----------------
    print("Writing output files...")

    # Sort samples by index (preserves samplelist order if provided)
    sorted_samples = sorted(allsampleindex.keys(), key=lambda k: allsampleindex[k])

    for chrom, regions in Rmap.items():
        output_filename = f"{args.output}.{chrom}.txt"
        # Sort regions by (start, end)
        Rlist = sorted(regions.items(), key=lambda item: (item[1][0], item[1][1]))

        with open(output_filename, "wt") as out:
            # Header line, compatible with DistractTree
            out.write(f"#{ref_label}\tChrom\tStart\tEnd\t" + "\t".join(sorted_samples) + "\n")

            for Rname, data in Rlist:
                start_trim = data[0]
                end_trim   = data[1]
                qnames_list = data[2]

                if not qnames_list:
                    continue

                Rinfo = [chrom, str(start_trim), str(end_trim)]

                # Initialize all sample cells as "NA"
                outList = ["NA"] * len(allsampleindex)

                # Group qnames by sample for this region
                sample_to_qnames = defaultdict(list)
                for q in qnames_list:
                    s = get_sample_name(q)
                    sample_to_qnames[s].append(q)

                # Decide cell content per sample
                for s_name, idx in allsampleindex.items():
                    q_list = sample_to_qnames.get(s_name)
                    if not q_list:
                        continue  # stays NA

                    # Exclude all firstcol alleles
                    candidates = [q for q in q_list if q not in excludelist]
                    if not candidates:
                        # All available alleles are excluded -> NA
                        continue

                    base = candidates[0]

                    # Attach any mapped alleles to this base if present in Converlist
                    if base in Converlist:
                        attached = Converlist[base]
                        cell_value = base + "," + ",".join(attached)
                    else:
                        cell_value = base

                    outList[idx] = cell_value

                out.write(Rname + "\t" + "\t".join(Rinfo + outList) + "\n")

    print("Done.")


if __name__ == "__main__":
    main()

