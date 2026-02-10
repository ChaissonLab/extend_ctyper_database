#!/usr/bin/env python3

import os
import argparse
import collections as cl
import numpy as np

nametochr = {'NC_060925.1': 'chr1', 'NC_060926.1': 'chr2', 'NC_060927.1': 'chr3', 'NC_060928.1': 'chr4', 'NC_060929.1': 'chr5', 'NC_060930.1': 'chr6', 'NC_060931.1': 'chr7', 'NC_060932.1': 'chr8', 'NC_060933.1': 'chr9', 'NC_060934.1': 'chr10', 'NC_060935.1': 'chr11', 'NC_060936.1': 'chr12', 'NC_060937.1': 'chr13', 'NC_060938.1': 'chr14', 'NC_060939.1': 'chr15', 'NC_060940.1': 'chr16', 'NC_060941.1': 'chr17', 'NC_060942.1': 'chr18', 'NC_060943.1': 'chr19', 'NC_060944.1': 'chr20', 'NC_060945.1': 'chr21', 'NC_060946.1': 'chr22', 'NC_060947.1': 'chrX', 'NC_060948.1': 'chrY'}


def main(args):
	
	

	with open(args.ref, mode = 'r') as f:
		
		reads = [x.splitlines() for x in f.read().split(">")[1:]]
		
		reads = {read[0].split()[0]: "".join(read[1:]) for read in reads}
		
	newfiles = set()
	index = 0
	with open(args.input, mode = 'r') as f:
		
		for line_ in f:
		
			index += 1	
			line = line_.strip().split()
			chr, start, end, blockname, name = line[0], int(line[1]), int(line[2]), line[-2] , line[-1]

			filename = args.output+"/" +name+ "/"+ name+".fa"
			themode = 'a'
			if filename not in newfiles:
				themode = 'w'
				newfiles.add(filename)	

			#name = nametochr[chr]+"_"+str(index) 
			os.makedirs(args.output+"/" +name, exist_ok=True)
			title = ">"+blockname + "\t" + line_.strip()
			
			with open(filename, mode = themode) as f:
				
				f.write(title+"\n")
				f.write(reads[chr][(start):(end)]+"\n")
				
			



def run():
	"""
		Parse arguments and run
	"""
	parser = argparse.ArgumentParser(description="program distract overlap genes and alignments on contigs")
	parser.add_argument("-i", "--input", help="path to input data file", dest="input", type=str, required=True)
	parser.add_argument("-r", "--ref", help="path to ref file", dest="ref", type=str, required=True)
	parser.add_argument("-o", "--out", help="path to output folder", dest="output", type=str, required=True)
	parser.set_defaults(func=main)
	args = parser.parse_args()
	args.func(args)
	
	
if __name__ == "__main__":
	run()
