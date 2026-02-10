#!/usr/bin/env python3

import argparse
import collections as cl
import gzip


def overlaploci(regions):
    """
    Given a list of regions: [ [start, end, ...], ... ],
    return a list of non-overlapping segments:
        [sub_start, sub_end, region_index]
    where each segment is covered by the region with the
    largest length among all overlapping regions.
    """
    if not regions:
        return []

    events = []
    sizes = []
    for idx, (start, end, *_) in enumerate(regions):
        events.append((start, 0, idx))  # 0 = start
        events.append((end,   1, idx))  # 1 = end
        sizes.append(end - start)

    # Sort by coordinate, then process starts before ends at same coord
    events.sort(key=lambda e: (e[0], e[1]))

    active = set()
    best_idx = None
    result = []

    for coord, etype, idx in events:
        if etype == 0:
            # region starts
            active.add(idx)
            # if this region is longer than current best, switch
            if best_idx is None or sizes[idx] > sizes[best_idx]:
                # close previous best segment at this coord
                if best_idx is not None:
                    result[-1][1] = coord
                # open new segment for this best region
                result.append([coord, regions[idx][1], idx])
                best_idx = idx

        else:
            # region ends
            if idx in active:
                active.remove(idx)

            if best_idx == idx:
                # best region ended, choose a new best if any
                if active:
                    new_best = max(active, key=lambda j: sizes[j])
                    if new_best != best_idx:
                        # close current segment at coord
                        result[-1][1] = coord
                        # open new best segment
                        result.append([coord, regions[new_best][1], new_best])
                        best_idx = new_best
                    else:
                        # same best, extend to coord
                        result[-1][1] = coord
                else:
                    # no regions left active, just close current segment
                    result[-1][1] = coord
                    best_idx = None

    return result


def open_maybe_gzip(path):
    """Open plain text or .gz transparently."""
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def regionindex(inputfile, output):
    """
    Parse headers like:
      >NAME1\tchr1:100-200+;NAME2\tchr1:150-250-;...
    Group by sample+contig, collapse overlaps with overlaploci(),
    and write:
      name  contig  orig_start  orig_end  sub_start  sub_end
    """
    # allregions[sample][contig] -> list of [start, end, name]
    allregions = cl.defaultdict(lambda: cl.defaultdict(list))

    with open_maybe_gzip(inputfile) as f:
        for line in f:
            if not line.startswith(">"):
                continue

            # Multiple regions may be on one header line separated by ';'
            theregions = line.strip().split(";")
            for theregion in theregions:
                if not theregion:
                    continue

                parts = theregion.split("\t")
                if len(parts) < 2:
                    continue

                name, region = parts[:2]

                strd = '+'
                # Drop trailing + / - if present
                if region and region[-1] in "+-":
                    strd = region[-1]
                    region = region[:-1]

                if ":" not in region or "-" not in region:
                    continue  # malformed region, skip

                contig, coord_str = region.split(":", 1)
                start_str, end_str = coord_str.split("-", 1)

                try:
                    start = int(start_str)
                    end = int(end_str)
                except ValueError:
                    continue

                # Define "sample" key from name: join fields 1 and 2
                # e.g. ">" + sample + "_" + something + ...
                name_clean = name.lstrip(">")
                name_parts = name_clean.split("_")
                if len(name_parts) >= 2:
                    key_name = "_".join(name_parts[1:3])
                else:
                    key_name = name_clean

                allregions[key_name][contig].append((start, end, name_clean+strd))

    with open(output, "w") as w:
        for sample, sample_regions in allregions.items():
            for contig, regions in sample_regions.items():
                newregions = overlaploci(regions)
                
                newregions_sort = [[] for x in regions]
                for region in newregions:
                    newregions_sort[region[-1]] = region
                
                for x,y in zip(newregions_sort, regions):
                    if len(x) == 0:
                        sub_start, sub_end = -1, -1
                    else:
                        sub_start, sub_end, idx = x
                    start, end, name_ = y
                    strd,name = name_[-1],name_[:-1]
                    # name, contig, original start, original end, sub_start, sub_end
                    w.write(
                        f"{contig}\t{start}\t{end}\t{name}\t{sub_end - sub_start}\t{strd}\t{sub_start - start if sub_start != -1 else sub_start}\t{end - sub_end if sub_end != -1 else sub_end}\n"
                    )

def main(args):
    regionindex(args.input, args.output)


def run():
    parser = argparse.ArgumentParser(
        description="Index and collapse overlapping regions from FASTA headers."
    )
    parser.add_argument(
        "-i", "--input",
        dest="input",
        required=True,
        help="Input file (plain text or .gz) with '>' region lines"
    )
    parser.add_argument(
        "-o", "--output",
        dest="output",
        required=True,
        help="Output TSV file"
    )

    parser.set_defaults(func=main)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    run()
