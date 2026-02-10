#!/usr/bin/env python3


import os
import collections as cl
import argparse
import re



def stretchresults(inputfile):
	
	with open(inputfile, mode ='r') as f:
		read = f.read()
		
	if read=='':
		return [0,"","",""]
	
	alllines=[a for a in read.split('\n')]
	#ref: original genome, query:assemblies
	alllinesindex=[a for a in range(len(alllines)) if len(alllines[a])>0 and alllines[a][0]!='#' and alllines[a][0]!=':' and  re.match(r'\S+\s\S+',alllines[a])!=None]
	
	eachline=alllinesindex[::2]
	qlines=[a for a in eachline]
	rlines=[a+2 for a in eachline]
	alines=[a+1 for a in eachline]
	if qlines==[]:
		return [-1,-1,-1,""]
	
	line0=alllines[qlines[0]]
	lineend = len(line0.strip())
	linestart = len(line0.split()[0]) + 1
	
	query=''.join([alllines[a][linestart:lineend].strip() for a in qlines])
	ref=''.join([alllines[a][linestart:lineend].strip() for a in rlines])
	align=''.join([(alllines[a])[linestart:lineend] for a in alines])[:len(ref)]
	
	
	return query, ref, align

def maskfile(inputfile, readsfile = ""):
	
	
	if len(readsfile) == 0:
		with open(inputfile, mode = 'r') as f:
			align = f.read()
		query = re.search(r'(?<=\[-asequence\]\s).+', align).group().strip()
		ref = re.search(r'(?<=\[-bsequence\]\s).+', align).group().strip()
	readsfile = query + "," + ref
	
	with open(readsfile.split(",")[0], mode = 'r') as f:
		
		query = [0,"".join(f.read().splitlines()[1:])]
		
	with open(readsfile.split(",")[1], mode = 'r') as f:
		ref = [0,"".join(f.read().splitlines()[1:])]

	
	output = []
	with open(inputfile, mode = 'r') as f:
		
		lineindex = 0
		for line in f:
			
			if len(line.strip()) == 0 or line.startswith(" ") or line.startswith('\t') or line.strip().startswith("#"):

				output.append(line)
				continue

			lineindex += 1
			
			if lineindex % 2 == 1:
				currentseq = query
			else:
				currentseq = ref

			consective= 0
			spacecount = 0
			for i,char in enumerate(line):
				if char not in [' ','\t']:
					if not consective:
						spacecount+=1
						if spacecount == 2:
							break				
					consective = 1
				else:
					consective = 0

			newline = line[:i]
			for char in line[i:]:
				
				if char.isalpha():
					newline += currentseq[1][currentseq[0]]
					currentseq[0] +=1

				else:
					newline += char
					
			output.append(newline)
		
	
	with open(inputfile, mode = 'w') as f:
		
		f.write("".join(output)+"\n")

def main(args):
	
	maskfile(args.input, args.reads)
		
def run():
	"""
		Parse arguments and run
	"""
	parser = argparse.ArgumentParser(description="program determine psuedogene")
	parser.add_argument("-i", "--input", help="path to input data file", dest="input", type=str, required=True)
	parser.add_argument("-r", "--reads", help="path to input data file", dest="reads", type=str, default="")
	parser.set_defaults(func=main)
	args = parser.parse_args()
	args.func(args)
	
	
if __name__ == "__main__":
	run()
