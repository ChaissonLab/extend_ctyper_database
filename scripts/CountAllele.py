#!/usr/bin/env python3

import collections as cl
import os
import argparse
from multiprocessing import Pool, cpu_count

# Shared globals for read-only use
name_to_allele = None
allele_groups = None

def loadtable(tablefile):
    allele_groups = cl.defaultdict(list)
    name_to_allele = {}

    with open(tablefile, mode='r') as f:
        for line in f:
            line = line.strip().split()
            if not line:
                continue
            alignto_group = line[-1].split("_")[0]
            alignto = alignto_group+"_"+line[0]
            names = line[-1].split(",")
            allele_groups[alignto_group].append(alignto)

            for name in names:
                name_to_allele[name] = alignto

    return name_to_allele, allele_groups


def init_globals(local_name_to_allele, local_allele_groups):
    global name_to_allele, allele_groups
    name_to_allele = local_name_to_allele
    allele_groups = local_allele_groups


def countsample(inputfile):
    global name_to_allele, allele_groups
    outputfile = inputfile + "_allele.out"
    allele_counts = cl.defaultdict(int)
    groupname = ""

    nonerrorgroup = set()
    errorgroup = set()
    with open(inputfile, mode='r') as f:
        for line in f:
            if len(line.strip()) == 0 or ("result:" not in line and line[0] != '#'):
                continue

            if line[0] == '>':
                groupname = line.strip().split()[0][2:]
                continue

            line = line.split("result:")[1]
            genes = [x for x in line.strip().split(",")[:-1]]

            if len(genes) == 0:
                errorgroup.add(groupname)
            else:
                nonerrorgroup.add(groupname)
                types = [name_to_allele.get(x, "NA") for x in genes]
                for thetype in types:
                    if thetype != "NA":
                        allele_counts[thetype] = types.count(thetype)


    for groupname in errorgroup - nonerrorgroup:
        alleles = allele_groups.get(groupname, [])
        for allele in alleles:
            allele_counts[allele] = 0
        


    with open(outputfile, mode='w') as f:
        for allele, count in allele_counts.items():
            if count == 0:
                f.write("{}\t-1\n".format(allele))
            else:
                f.write("{}\t{}\n".format(allele, count))


def main(args):
    local_name_to_allele, local_allele_groups = loadtable(args.table)

    if os.path.isdir(args.input):
        allfiles = [os.path.join(args.input, x) for x in os.listdir(args.input) if not x.endswith("_allele.out") and x + "_allele.out" not in os.listdir(args.input)]
    else:
        allfiles = [args.input]

    with Pool(processes=args.threads, initializer=init_globals,
              initargs=(local_name_to_allele, local_allele_groups)) as pool:
        pool.map(countsample, allfiles)


def run():
    parser = argparse.ArgumentParser(description="Determine allele-specific gene counts from ctyper result")
    parser.add_argument("-i", "--input", help="Path to input data file or directory", dest="input", type=str, required=True)
    parser.add_argument("-t", "--table", help="Path to annotation table file", dest="table", type=str, required=True)
    parser.add_argument("-o", "--output", help="(Unused, kept for compatibility)", dest="output", type=str, default="")
    parser.add_argument("-n", "--threads", help="Number of threads (default: all cores)", dest="threads", type=int, default=cpu_count())
    parser.set_defaults(func=main)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    run()

