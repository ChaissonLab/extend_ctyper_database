#!/usr/bin/env python3

import os
import re
import collections as cl
import math
import numpy as np
import random as rd
import argparse
import gzip

class Matrix:
	
	def __init__(self):
		
		self.dim = 0
		self.size = 0
		self.data = np.zeros((0,0), dtype=np.int32)
		
	def load(self, filepath):
		
		with gzip.open(filepath,'rb') as f:
			
			firstrow = f.readline().decode().strip().split(",")
			if firstrow[-1] == "":
				firstrow = firstrow[:-1]
				
			firstrow = [int(x) for x in firstrow]
			
			self.dim = len(firstrow)
			
			self.data = np.zeros((self.dim,self.dim), dtype=np.int32)
			
			self.data[0,:] = firstrow 
			
			index = 0
			for row in f:
				
				index += 1
				row = row.decode().strip().split(",")[:-1]
				self.data[index,:] = [0]*(self.dim- len(row))+[int(x) for x in row]
				
			for i in range(self.dim):
				
				for j in range(i+1, self.dim):
					
					self.data[j,i] = self.data[i,j]
					
		return self
	
def kmerencode(kmer):
	
	kmerint = 0
	for base in kmer:
		
		kmerint *= 4
		if base=='a' or base =='A':
			
			kmerint += 0
			
		elif base=='c' or base =='C':
			
			kmerint += 1
			
		elif base=='g' or base =='G':
			
			kmerint += 2
		elif base=='t' or base =='T':
			
			kmerint += 3
			
	return kmerint

def kmerdecode(kmerint, size = 31):
	
	index = 0
	text = ['A' for a in range(size)]
	
	while kmerint:
		
		text[index] = "ACGT"[kmerint % 4]
		kmerint //= 4
		index += 1
		
	return "".join(text[::-1])


def intencode(value):
	
	code = ""
	
	code = chr( ord('0') + (value % 64) ) + code
	value = value//64
	
	while value:
		
		code = chr( ord('0') + (value % 64) ) + code
		value = value//64
		
	return code


def find_identical_rows(matrix: np.ndarray):
	row_map = {}
	unique_rows = []
	index_groups = {}
	exclude_indeces = []
	
	for i, row in enumerate(matrix):
		row_tuple = tuple(np.round(row, decimals=8))  # Use rounding to avoid float noise
		if row_tuple in row_map:
			row_map[row_tuple].append(i)
		else:
			row_map[row_tuple] = [i]
			
	for indices in row_map.values():
		exclude_indeces.extend(indices[1:])  # Take the first as representative
		index_groups[indices[0]] = indices
		
	return exclude_indeces, index_groups

class KmerData:
	
	def __init__(self, seqfile,  kmersize = 31, msize = 0, leave = []):
		
		
		
		self.exclude = set(leave)
		
		self.kmersize = kmersize
		
		self.kmerslist = dict()
		
		self.genekmercounts = []
		
		self.samplesizes = []   
		
		index = 0
		name = ""
		samples = []
		with open(seqfile,mode = 'r') as r:
			for line in r:
				
				line = line.strip()
				
				if len(line):
					if line[0]==">":
						name = line.split()[0][1:]
						if len([x for x in self.exclude if len(x) and x in name]):
							name = ""
							continue
						self.samplesizes.append(0)
						samples.append(name)
						self.genekmercounts.append( 0 )
						index += 1
					else:
						if name != "":
							self.samplesizes[-1] += len(line.strip())
							
							
		self.sampleslist = list(samples)
		
		
		if msize < 20000000:
			self.matrix = kmer_cluster32.SparseKmerMartrix(len(self.sampleslist), 0)
		else:
			self.matrix = kmer_cluster32.SparseKmerMartrix(len(self.sampleslist), 0)
			
	def LoadKmers(self, kmerfile):
		
		kmerindex = 0
		with open(kmerfile, mode = 'r') as f:
			
			for line in f:
				
				if line.startswith(">"):
					
					self.matrix.addkmer()
					self.kmerslist[kmerencode(line.strip().split()[0])] = kmerindex
					kmerindex += 1
					
					
					
	def AddKmer(self,kmer, index):
		
		kmerindex = self.kmerslist.get(kmer, -1) 
		
		if kmerindex >= 0:
			
			self.genekmercounts[index] += 1
			
			self.matrix.add(kmerindex, index)
			
		return kmerindex
	
	
	def ReadKmers(self, seqfile, outputfile):
		
		intoperator = 4**(self.kmersize-1)
		
		current_k = 0
		reverse_k = 0
		current_size = 0
		
		qname = ""
		
		currentline = ""
		
		with open(outputfile,mode = 'w') as w:
			
			with open(seqfile,mode = 'r') as r:
				
				posi_linestart = 0
				index = -1
				for line in r:
					
					if len(line) ==0:
						continue
					
					if line[0]==">":
						current_k=0
						reverse_k=0
						current_size = 0
						posi_linestart = 0
						
						qname = line.strip().split()[0][1:]
						if len([x for x in self.exclude if len(x) and x in qname]):
							qname = ""
							continue
						
						index += 1
						
						w.write(line)
						continue
					
					if qname == "":
						continue
					
					line_out = bytearray(line.lower().encode())
					for posi,char in enumerate(line.upper()):
						
						if char =="\n":
							continue
						
						if char not in ['A','T','C','G']:
							
							current_k = 0
							reverse_k = 0
							current_size = 0
							
							continue
						
						if current_size >= self.kmersize:
							
							current_k %= intoperator
							current_k <<= 2
							current_k += ['A','C','G','T'].index(char)
							
							reverse_k >>= 2
							reverse_k += (3-['A','C','G','T'].index(char)) * intoperator
							
						else:
							
							
							current_k <<= 2
							current_k += ['A','C','G','T'].index(char)
							reverse_k += (3-['A','C','G','T'].index(char)) * ( 1 << (2*current_size))
							
							current_size += 1
							
						if current_size >= self.kmersize :  
							
							kmer = max(current_k, reverse_k)
							
							findindex = self.AddKmer(kmer,index)
							
							if findindex >= 0:
								
								line_out[posi] = line_out[posi] - 32
								
					posi_linestart += len(line)
					
					line_out = line_out.decode()
					
					w.write(line_out)
					
					
def sample1000(thelist):
	
	if len(thelist) < 1000:
		return thelist
	else:
		ratio = 1000/len(thelist)
		thelist = [x for i,x in enumerate(thelist) if int(i * ratio + 0.5)!=int((i-1) * ratio + 0.5)]
		return thelist
	
def build_balanced_newick(leaves, up_dist=0.0) -> str:
	"""
	Recursively builds a balanced binary tree from a list of leaves.
	Assigns 0 branch length to inner clades, and total `up_dist` to the root.
	"""
	if len(leaves) == 1:
		return leaves[0] + ":" + '{:.7f}'.format(up_dist)
	elif len(leaves) == 2:
		return f"({leaves[0]}:0.0000000,{leaves[1]}:0.0000000):{up_dist:.7f}"
	else:
		mid = len(leaves) // 2
		left = build_balanced_newick(leaves[:mid], 0.0)
		right = build_balanced_newick(leaves[mid:], 0.0)
		return f"({left},{right}):{up_dist:.7f}"
	
class UPGMANode:
	
	def __init__(self, left=None, right=None, up_dist=0.0, down_dist=0.0, index = None):
		
		self.index = index
		self.left = left
		self.right = right
		self.up_dist = up_dist
		self.down_dist = down_dist
		self.nleaves = 1
		
		if type(left) is type(self) and type(right) is type(self):
			self.nleaves += left.nleaves + right.nleaves
		else:
			self.nleaves = 1
			
			
	def leaves(self) -> list:
		
		if self is None:
			return []
		elif self.right == None:
			return [self.left]
		
		return self.left.leaves() + self.right.leaves() 
	
	def to_newick(self, group_dict = dict()) -> str:
		
		if self.right == None:
			
			name = self.left
			
			if name in group_dict and len(group_dict[name]) > 1:
				
				return build_balanced_newick(group_dict[name], self.up_dist)
			else:
				return name + ":" + '{:.7f}'.format(self.up_dist)
			
			return self.left + ":" + '{:.7f}'.format(self.up_dist)
		else:
			return ( "(" + ",".join([x.to_newick(group_dict) for x in [self.left, self.right]]) + "):"+ '{:.6f}'.format(self.up_dist) )
		
	def allleaves(self):
		
		if type(self.left) is not type(self) or  type(self.right) is not type(self):
			return [self]
		
		
		return self.left.allleaves() + self.right.allleaves()
	
	def reorder(self, dist_matrix, sibling_indexes = None, sign = 0) -> list:
		
		if type(self.left) is not type(self) or  type(self.right) is not type(self):
			return
		
		leftleaves_indexes = [x.index for x in self.left.allleaves()]
		rightleaves_indexes = [x.index for x in self.right.allleaves()]
		
		self.left.reorder(dist_matrix, rightleaves_indexes, 1)
		self.right.reorder(dist_matrix, leftleaves_indexes, -1)
		
		if sibling_indexes == None or sign == 0:
			return 
		
		left_distance_mean = np.mean([ dist_matrix[i,j] for j in sample1000(sibling_indexes) for i in sample1000(leftleaves_indexes)])
		
		right_distance_mean = np.mean([ dist_matrix[i,j] for j in sample1000(sibling_indexes) for i in sample1000(rightleaves_indexes)])
		
		if sign * right_distance_mean > sign * left_distance_mean:
			
			temp = self.left
			self.left = self.right
			self.right = temp
			
			
			
class UPGMA:
	
	
	def __init__(self, dist_matrix: np.ndarray, header: list):
		
		self.distances = dist_matrix
		self.header = header
		self.exclude = [0 for x in range(len(header))]+[1 for x in range(len(header))]
		self.build_tree(self.distances, self.header)
		
	def getmindist(self)-> tuple:
		
		min_row, min_col = 0,0
		
		curr_min = np.inf
		for row,values in enumerate(self.work_matrix):
			
			if self.exclude[row]:
				continue
			
			for col, value in enumerate(values):
				
				if self.exclude[col]:
					continue
				
				if value == 0:
					
					return (row, col)
				
				elif value < curr_min:
					
					curr_min = value
					min_row = row
					min_col = col
					
					
		return (min_row, min_col)
	
	
	def build_tree(self, dist_matrix: np.ndarray, header: list) -> UPGMANode:
		
		nodes = [UPGMANode(taxon, index = i) for i, taxon in enumerate(header)] 
		
		headertoindex = {name:i for i,name in enumerate(header)}
		
		self.size = 2*len(header)
		
		self.work_matrix = np.array([row+[np.inf]*len(row) for row in dist_matrix.tolist()]+[[np.inf]*(2*len(row)) for row in dist_matrix.tolist()], dtype=float)
		np.fill_diagonal(self.work_matrix, np.inf)
		
		new_node = None 
		for turn in range(len(nodes)-1):
			
			least_id = self.getmindist()
			least_dist = self.work_matrix[least_id[0], least_id[1]]
			node1, node2 = nodes[least_id[0]], nodes[least_id[1]]
			self.exclude[least_id[0]] = 1
			self.exclude[least_id[1]] = 1
			
			new_node = UPGMANode(node2, node1)
			nodes.append(new_node)
			self.exclude[len(nodes)-1] = 0
			node1.up_dist = least_dist / 2 - node1.down_dist
			node2.up_dist = least_dist / 2 - node2.down_dist
			new_node.down_dist = least_dist / 2
			
			# create new working distance matrix
			self.update_distance( nodes, least_id)
			
		self.tree = new_node 
		
	def update_distance(self, nodes: list, least_id: tuple) -> np.ndarray:
		
		length = len(nodes)
		nleaves1, nleaves2 = nodes[least_id[0]].nleaves, nodes[least_id[1]].nleaves
		nleaves  = nleaves1 + nleaves2 
		
		for i in range(length-1):
			
			if self.exclude[i]:
				continue
			
			
			self.work_matrix[i,length-1] = ( self.work_matrix[i][least_id[0]]*nleaves1 + self.work_matrix[i][least_id[1]]*nleaves2 ) / nleaves
			self.work_matrix[length-1,i] = self.work_matrix[i,length-1]
			
def projectiontree(matrix, header, exclude = set(), index_groups = dict()):
	
	if len(header) < 2:
		return header[0] + ";"
	
	matrix = np.matrix(matrix)
	
	index_groups_header = {header[k]:[header[x] for x in v] for k,v in index_groups.items()}
	header_used = [x for i,x in enumerate(header) if i not in exclude]
	useindex= [i for i,x in enumerate(header) if i not in exclude]
	
	dm = np.ones((len(useindex),len(useindex)))
	vars = [math.sqrt(matrix[i,i]) for i in range(len(matrix))]
	for i,x in enumerate(useindex):
				
		var1 = max(1.0,vars[x])
		
		for j,y in enumerate(useindex[:i]):
			
			var2 = max(1.0,vars[y])
			dm [i,j] = 1 - matrix[x,y]/(var1 *var2)
			dm [j,i] =  dm [i,j]
			
	tree = UPGMA(dm, header_used).tree
	tree.reorder(dm)
	
	return tree.to_newick(index_groups_header)+";"


def ReorderFile(infile,outfile,treeorder, name = ""):
	
	with open(infile,mode = 'r') as r:
		
		reads = r.read().split(">")[1:]
		
	reordered = [reads[index].splitlines() for index in treeorder]
	
	if len(name):
		reordered = [name+read[0].replace(" ","\t")+"\n"+"\n".join(read[1:]) for read in reordered]
	else:
		reordered = [read[0].replace(" ","\t")+"\n"+"\n".join(read[1:]) for read in reordered]
		
	with open(outfile,mode = 'w') as w:
		
		w.write(">"+"\n>".join(reordered)+"\n")
		
	header = [read.splitlines()[0].replace(" ","\t") for read in reordered]
	
	return header

def ReorderMatrix(matrix, order):
	
	
	matrix_rshape = np.reshape(matrix, (len(order), len(order)))
	
	#order_sort = sorted(range(len(order)), key = lambda x: order[x])
	
	newmatrix = []
	for reorderindex, index in enumerate(order):
		
		row = matrix_rshape[index].tolist()
		row_rshape = [row[i] for i in order]
		
		newmatrix.append(row_rshape)
		
	return newmatrix


def MatrixFulltoUpper(matrix, size):
	
	newmatrix = []
	
	for i in range(size):
		
		newmatrix += matrix[(i*size + i) : (i*size + size) ] 
		
	return newmatrix

def ReorderNewickname(text):
	
	locations = [x.span() for x in re.finditer(r'[^(^)^:^,^;]+', text) if text[x.span()[0]-1] != ":"]
	
	last_location = 0
	newtext = ""
	nameindex = 0
	for location in locations:
		
		newtext += text[last_location:location[0]] + str(nameindex)+"_"
		
		last_location = location[1]
		
		nameindex += 1
		
	newtext += text[last_location:]
	
	return newtext

def hashrow(kmerindex , genenum, KmerReader ):
	
	row = tuple(sorted(list(KmerReader.matrix.getkmerrow( kmerindex ) ) ) )
	
	return ( - (abs( len(row) -  genenum/2 - 0.1 ))  , hash(row) )


def getcontigs(inputfile):
	
	names = []
	
	with open(inputfile, mode = 'r') as f:
		
		for line in f:
			
			if line.startswith(">"):
				name = line[1:].split()[0]
				names.append(name)
				
	return names



def annotatefasta(seqfile, normfile, kmerfile, outputfile, leave = []):
	
	if len(normfile) == 0:
		
		import pyximport; pyximport.install(setup_args={"include_dirs":np.get_include()})
		import kmer_cluster
		import kmer_cluster32
		
		kfilesize = os.stat(kmer).st_size
		
		KmerReader = KmerData(seqfile, msize = kfilesize, leave = leave)
		
		KmerReader.sampleslist = [name+x for x in KmerReader.sampleslist]
		samplelist = KmerReader.sampleslist
		
		KmerReader.LoadKmers(kmerfile)
		KmerReader.ReadKmers(seqfile,outputfile)
		
		sqmatrix = np.reshape( KmerReader.matrix.SquareMatrix(0), (len(samplelist), len(samplelist)) )
	else:
		
		sqmatrix = Matrix().load(normfile).data
		samplelist = getcontigs(seqfile)
	
	
	exclude_indeces, index_groups = find_identical_rows(sqmatrix)
		
	treetext = projectiontree( sqmatrix, samplelist, exclude_indeces, index_groups)
	
	treeorder = [x.split(":")[0].replace("(", "") for x in treetext[:-1].split(",")]

	nametoindex = {x:i for i,x in enumerate(samplelist)}
	
	treeorder = [nametoindex[sample] for sample in treeorder]
	
	headers = ReorderFile(seqfile, outputfile, treeorder)
	
	newnorm = np.array(ReorderMatrix(sqmatrix, treeorder))
	
	return headers,treetext,newnorm 


def main(args):
	
	if len(args.output) == 0:
		args.output = args.seq + "_treeorder.fa"
		
	headers, treetext, newnorm = annotatefasta(args.seq, args.norm, args.kmer, args.output, args.leave.split(","))
	
	with open(args.output + "_tree.ph", mode = 'w') as f:
		
		f.write(treetext + "\n")
		
	with gzip.open(args.output + "_norm.txt.gz", "wt") as f:
		np.savetxt(f, newnorm, delimiter=',', fmt='%1.4f')
		
def run():
	"""
		Parse arguments and run
	"""
	parser = argparse.ArgumentParser(description="program determine tree")
	parser.add_argument("-i", "--input", help="path to input data file",dest="seq", type=str, required = True)
	parser.add_argument("-n", "--norm", help="path to input data file",dest="norm", type=str, default = '')
	parser.add_argument("-k", "--kmer", help="path to input data file",dest="kmer", type=str, default = '')
	parser.add_argument("-o", "--output", help="path to input data file", dest="output", type=str, default = "")
	parser.add_argument("-l", "--leave", help="path to input data file", dest="leave", type=str, default = '')
	parser.set_defaults(func=main)
	args = parser.parse_args()
	args.func(args)
	
	
if __name__ == "__main__":
	run()
