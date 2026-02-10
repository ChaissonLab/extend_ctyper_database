#!/usr/bin/env python3

import os
import argparse
import collections as cl
import numpy as np

nametochr = {'NC_060925.1': 'chr1', 'NC_060926.1': 'chr2', 'NC_060927.1': 'chr3', 'NC_060928.1': 'chr4', 'NC_060929.1': 'chr5', 'NC_060930.1': 'chr6', 'NC_060931.1': 'chr7', 'NC_060932.1': 'chr8', 'NC_060933.1': 'chr9', 'NC_060934.1': 'chr10', 'NC_060935.1': 'chr11', 'NC_060936.1': 'chr12', 'NC_060937.1': 'chr13', 'NC_060938.1': 'chr14', 'NC_060939.1': 'chr15', 'NC_060940.1': 'chr16', 'NC_060941.1': 'chr17', 'NC_060942.1': 'chr18', 'NC_060943.1': 'chr19', 'NC_060944.1': 'chr20', 'NC_060945.1': 'chr21', 'NC_060946.1': 'chr22', 'NC_060947.1': 'chrX', 'NC_060948.1': 'chrY'}

def mergeregions(coordis):
	
	names = [y[2] for y in coordis]
	coordis = [x for y in coordis for x in y[:2]]
	
	coordis_sort = sorted(range(len(coordis)), key = lambda x:  coordis[x])
	
	genes = set()
	segments = []
	cover = 0
	start = 0
	for i,x in enumerate(coordis_sort):
		
		coordi = coordis[x]
		name = names[x//2]
		
		genes.add(name)
		
		if x % 2:
			cover -= 1
			if cover == 0:
				segments.append([start, coordi, ";".join(genes)])
				genes = set()
		else:
			cover += 1
			if cover == 1:
				start = coordi
	
	
	return segments

def readgff3(bedfile):
	
	genes_region = cl.defaultdict(list)
	with open(bedfile, mode = 'r') as f:
		
		for line in f:
			
			if line.startswith('#'):
				continue
			
			line = line.strip().split()
			
			
			if line[2] != "exon" :
				continue
			
			info = line[-1].split(";")
			name = [x for x in info if x.startswith("source_gene=")][0][12:]
			
			name = "-".join(name.split("-")[:1]) if "-" in name else name
			
			start = int(line[3])
			end = int(line[4])
			
			genes_region[line[0]].append([start,end, name])

	genes_region = {chr:mergeregions(region) for chr,region in genes_region.items()} 
			
	return genes_region

def main(args):
	
	genes_region = readgff3(args.genecode)
	
    
	
	with open(args.ref, mode = 'r') as f:
		
		reads = [x.splitlines() for x in f.read().split(">")[1:]]
		
		reads = {read[0].split()[0]: "".join(read[1:]) for read in reads}
		
	newfiles = set()	
	index = 0
	with open(args.input, mode = 'r') as f:
		
		for line_ in f:
			index += 1
			line = line_.strip().split()
			chr, start, end,  blockname, name = line[0], int(line[1]), int(line[2]),  line[-2], line[-1]

			chr_ = nametochr[chr]
                
			filename = args.output+"/"+name + "/"+ name+"_exons.fa"
			themode = 'a'
			if filename not in newfiles:
				newfiles.add(filename)	
				themode = 'w'					

			genes_region_chr = genes_region[chr_]
			os.makedirs(args.output+"/" +name, exist_ok=True)
			gene_overlaps = [x for i,x in enumerate(genes_region_chr) if max(x[1], end) - min(x[0],start) - (end - start) - (x[1] - x[0]) < 0]
			index2 = 0	
			with open(filename, mode = themode) as f:
				
				for gene_overlap in gene_overlaps:
		
					index2 += 1	
					start0, end0, title = gene_overlap

					f.write(">Trans_{}_{}\t{}\t{}\t{}\t{}\n".format(blockname,index2,title, chr, max(start0,start), min(end0, end)))
					f.write(reads[chr][(max(start0,start)):(min(end0, end))]+"\n")
				
				
def run():
	"""
		Parse arguments and run
	"""
	parser = argparse.ArgumentParser(description="program distract overlap genes and alignments on contigs")
	parser.add_argument("-i", "--input", help="path to input data file", dest="input", type=str, required=True)
	parser.add_argument("-r", "--ref", help="path to ref file", dest="ref", type=str, required=True)
	parser.add_argument("-g", "--genecode", help="path to ref file", dest="genecode", type=str, required=True)
	parser.add_argument("-o", "--out", help="path to output folder", dest="output", type=str, required=True)
	parser.set_defaults(func=main)
	args = parser.parse_args()
	args.func(args)
	
	
if __name__ == "__main__":
	run()
