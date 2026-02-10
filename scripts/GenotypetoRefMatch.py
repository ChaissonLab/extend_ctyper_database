#!/usr/bin/env python3

import collections as cl

mainchr = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM']

if __name__ == "__main__":
    
    allele = "ALL_refmatch.txt"
    genotype = "HG002_leave.out"
    genotype = "oldctyper/ctypernocorr.out"
    #genotype = "NA18967.final.cram"
    genotype_out = "HG002.final.cram.genotype.refmatch2"
    assem_out = "HG002.final.cram.assem.refmatch2"
    hg38_out = "HG002.final.cram.ref.refmatch"
    pref = "HG002"
    
    groupcounts = "./CHM13num.txt"

    mainchrname = set()
    refmatch = cl.defaultdict(str)
    with open(allele, mode = 'r') as f:
        
        for line in f:
            line = line.split()
            refmatch[line[0]] = line[3]
            if line[2].split(":")[0] in mainchr:
                 mainchrname.add(line[0])


    groupnumbers = cl.defaultdict(int)
    with open(groupcounts, mode ='r') as f:
        for line in f:
            line = line.split()
            groupnumbers[line[0]] =int(line[1])

    redundants = set()
    passgroup = set()
    results = []
    with open(genotype, mode = 'r') as f, open(genotype_out, mode = 'w') as w:
        chm13counts = 0
        skip = 0
        for line in f:
            if line.startswith(">#"):
                groupname = line.split('\t')[0][1:]
                groupnum = groupnumbers[groupname]

                skip = 0
                if groupname in redundants:
                    skip = 1
                else:
                    redundants.add(groupname)
            
            if skip == 1:
                continue

            if line.startswith("result:") and line.count(",") >0 and ( line.count(",") < max(5*groupnum, groupnum + 5) and groupnum > 0 ):
                passgroup.add(groupname[1:])
                alleles = line.split()[1].split(",")[:-1]
                refmatches = ["{}\t{}".format(allele, refmatch[allele]) for allele in alleles]
                
                w.write("\n".join(refmatches)+"\n")

    with open(assem_out, mode = 'w') as w:
        for name, ref in refmatch.items():
            groupname = name.split("_")[0]
            if name.split("_")[1] == pref and groupname in passgroup :
                w.write("{}\t{}".format(name, ref)+"\n")

   
    with open(hg38_out, mode = 'w') as w:
        for name, ref in refmatch.items():
            groupname = name.split("_")[0]
            if name.split("_")[1] == "HG38" and groupname in passgroup and  name in mainchrname :
                w.write("{}\t{}".format(name, ref)+"\n")


