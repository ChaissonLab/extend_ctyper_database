#!/usr/bin/env python3


import re
import os
import collections as cl
import math
import numpy as np
import scipy 
import argparse
import gzip

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


def intencode(value, length = 0):
	
	code = ""
	
	code = chr( ord('0') + (value % 64) ) + code
	value = value//64
	
	while value:
		
		code = chr( ord('0') + (value % 64) ) + code
		value = value//64
	if len(code) < length :
		code = (length-len(code))*"0" + code
		
	return code
def MatrixFulltoUpper(matrix, size):
	newmatrix = []
	for i in range(size):
		newmatrix.append(matrix[(i*size + i) : (i*size + size) ] )
	return newmatrix
class UPGMANode:
	
	def __init__(self, left=None, right=None, up_dist=0.0, down_dist=0.0, index = None):
		
		self.index = index
		self.left = left
		self.right = right
		self.up_dist = up_dist
		self.down_dist = down_dist
		self.nleaves = 1
		
		if type(left) is type(UPGMANode):
			self.nleaves += left.nleaves + right.nleaves
		else:
			self.nleaves = 1
			
			
	def leaves(self) -> list:
		
		if self is None:
			return []
		elif self.right == None:
			return [self.left]
		
		return self.left.leaves() + self.right.leaves() 
	
	def to_newick(self) -> str:
		
		if self.right == None:
			return self.left + ":" + '{:.7f}'.format(self.up_dist)
		else:
			return ( "(" + ",".join([x.to_newick() for x in [self.left, self.right]]) + "):"+ '{:.7f}'.format(self.up_dist) )
		
	def allleaves(self):
		
		if type(self.left) is not type(UPGMANode) or  type(self.right) is not type(UPGMANode):
			return self
		
		
		return self.left.allleaves() + self.right.allleaves()
	
	def reorder(self, dist_matrix, sibling_indexes = None, sign = 0) -> list:
		
		if type(self.left) is not type(UPGMANode) or  type(self.right) is not type(UPGMANode):
			return
		
		leftleaves_indexes = [x.index for x in self.left.allleaves()]
		rightleaves_indexes = [x.index for x in self.right.allleaves()]
		
		self.left.reorder(dist_matrix, rightleaves_indexes, 1)
		self.right.reorder(dist_matrix, leftleaves_indexes, -1)
		
		if sibling_indexes == None or sign == 0:
			return 
		
		left_distances = [ dist_matrix[i,j] for j in sibling_indexes for i in leftleaves_indexes]
		left_distance_mean = sum(left_distances)/len(left_distances)
		
		right_distances = [ dist_matrix[i,j] for j in sibling_indexes for i in rightleaves_indexes]
		right_distance_mean = sum(right_distances)/len(right_distances)
		
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
			

def group_matrix(matrix, size, similar):

        groups = [[] for x in range(size)]
        group_edges = [x for x in range(size)]
        for i in range(size):

                indexstart = i*size

                diag_cutoff = matrix[indexstart + i ]* similar

                for j in range(i+1, size):

                        edge = matrix[indexstart + j ]


                        if edge > max ( diag_cutoff , matrix[j*size + j ]* similar ):

                                if group_edges[i] < group_edges[j]:
                                        group_edges[j] = group_edges[i]
                                else:   
                                        group_edges[i] = group_edges[j]

        for i, group_index in enumerate(group_edges):

                last_group_index = i

                while group_index < last_group_index:

                        last_group_index = group_index
                        group_index = group_edges[group_index]

                group_edges[i] = group_index
                groups[group_index].append(i)

        return groups, group_edges



def projectiontree(matrix, header):
	
	if len(header) < 2:
		return header[0] + ";"
	
	matrix = np.matrix(matrix)
	
	dm = np.ones((matrix.shape))
	vars = [matrix[i,i] for i in range(len(matrix))]
	for i in range(len(matrix)):
		
		var1 = max(1.0,vars[i])
		for j in range(i):
			
			var2 = max(1.0,vars[j])
			
			dm [i,j] = 1 - matrix[i,j]/math.sqrt((var1*var2))
			dm [j,i] =  dm [i,j]
			
	tree = UPGMA(dm, header).tree
	tree.reorder(dm)
	
	return tree.to_newick()+";"

def ReorderFile(seqfile, treeorder, name):
	
	with open(seqfile,mode = 'r') as r:
		
		reads = r.read().split(">")[1:]
		
	reordered = [reads[index].splitlines() for index in treeorder]
	reordered = [name+"_"+read[0].replace(" ","\t")+"\n"+"\n".join(read[1:]) for read in reordered]
	
	with open(seqfile,mode = 'w') as w:
		
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
		
		newmatrix.extend(row_rshape[reorderindex:])
		
	return newmatrix

def MatrixFilter(matrix, size, filterindex):
	
	filterindex = set(filterindex)
	
	newmatrix = []
	
	for i in range(size):
		
		if i not in filterindex:
			
			newmatrix.extend([matrix[i*size + j] for j in range(size) if j not in filterindex] )
			
	return newmatrix

def MatrixUppertoFull(matrix, size):
	newmatrix = []
	pos  = 0	
	for i in range(size):
		newmatrix.extend([ newmatrix[x*size + i] for x in range(i) ] )	
		newmatrix.extend(matrix[pos : (pos + size - i) ] )
		pos += size - i
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




class KmerData:
	
	def __init__(self, seqfile,  kmersize = 31, msize = 10):
		
		self.kmersize = kmersize
		
		self.kmerslist = dict()
		
		self.genekmercounts = []
		
		self.samplesizes = []  

		self.kmerflag = dict() 
			
		self.headers = []
		self.msize = msize	
		self.numpath = 0
		self.pathsize = cl.defaultdict(int)
		self.sampleloc = dict()	
	
		index = 0
		samples = []
		with open(seqfile,mode = 'r') as r:
			for line in r:
				
				line = line.strip()
				
				if len(line):
					if line[0]==">":
						self.headers.append(line)
						self.samplesizes.append(0)
						name = line.split()[0][1:]
						samples.append(name)
						self.sampleloc[name] = [name, line.split()[1]]
						if "_group" in name:
							self.sampleloc["_".join(name.split("_")[2:])] = [name, line.split()[1]]

						index += 1
					else:
						self.samplesizes[-1] += len(line.strip())
						
						
		self.sampleslist = list(samples)

		self.initiatematrix()		


	def initiatematrix(self):
		import pyximport; pyximport.install(setup_args={"include_dirs":np.get_include()})
		import kmer_cluster
		import kmer_cluster32
		self.genekmercounts = []

		if self.msize < 20000000 or 1==1:
			self.matrix = kmer_cluster32.SparseKmerMartrix(len(self.sampleslist), 0)
		else:
			self.matrix = kmer_cluster32.SparseKmerMartrix(len(self.sampleslist), 0)

		for i in range(len(self.kmerslist)):
			
			self.matrix.addkmer()
		
	
		
	def LoadKmers(self, kmerfile):
	
		headers = []	
		kmerindex = 0
		with open(kmerfile, mode = 'r') as f:
			
			for line in f:
				if line[0] == "+":
					if len(line):
						name = "_".join(line.split()[0].split('_')[:-2])[1:]
						used_region = [int(line.split()[0].split('_')[-2]), int(line.split()[0].split('_')[-1])]
						locus = self.sampleloc[name][1]
						contig = locus.split(":")[0]
						region = locus.split(":")[1][:-1]
						strd = locus[-1]
						start,end = int(region.split("-")[0]), int(region.split("-")[1])
						if used_region[1] == -1:
							used_region[1] = end - start

						if strd == '+':
							abso_region = [used_region[0] + start, used_region[1] + start]
						else:
							abso_region = [end - used_region[1], end - used_region[0]]

						header = "+{}\t{}-{}\t{}:{}-{}{}".format(self.sampleloc[name][0], used_region[0] , used_region[1] ,contig, abso_region[0], abso_region[1], strd)
					else:
						header = "+"+line[1:].strip()

					headers.append(header)
				
				elif len(line) and line[0] != ">":
					
					self.matrix.addkmer()
					
					elements = line.strip().split()
					
					if len(elements) > 2:
						kmer, pathindex, pathloc,flag,ratio = elements
						
						self.pathsize[pathindex] = max(self.pathsize[pathindex], int(pathloc) + 1)
						self.numpath = max(self.numpath, int(pathindex))
					else:
						kmer = elements[0]
						pathindex = "0"
						pathloc = "0"
						
					self.kmerslist[kmerencode(kmer)] = ( kmerindex , pathindex , pathloc, flag , ratio )
					kmerindex += 1
		
		return "\n".join(headers)
					
					
	def AddKmer(self,kmer, index):
		
		kmerindex = self.kmerslist.get(kmer, (-1,-1,-1,0) )[0] 
		
		if kmerindex >= 0:
	
			self.genekmercounts[index] += 1
			self.matrix.add(kmerindex, index)
			
		return kmerindex
	
	
	def ReadKmers(self, seqfile, filternames = set([])):
		
		intoperator = 4**(self.kmersize-1)
		
		current_k = 0
		reverse_k = 0
		current_size = 0
		
		currentline = ""
		
		iffilter = 0
		
		with open(seqfile,mode = 'r') as r:
			
			index = -1
			for line in r:
				
				if len(line) ==0:
					continue
				
				if line[0]==">":
					
					if line.split()[0][1:] in filternames:
						iffilter = 0
						continue
					else:
						iffilter = 1
						
						
					self.genekmercounts.append(0)
					index += 1
					current_k=0
					reverse_k=0
					current_size = 0
					
					posi_linestart = 0
					
					continue
				
				elif iffilter == 0:
					continue
				
				
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
						
						
class node:
	
	def __init__(self, parent = None, name = "", distance =0.0, index = 0):
		
		self.children = []
		
		self.parent = parent
		
		self.name = name
		
		self.distance = distance
		
		self.index = index
		
		self.annotation = 0
		
		if parent is not None:
			
			parent.children.append(self)
			
	def __str__(self):
		
		if len(self.children) == 0:
			
			return self.name+":"+"{:.7f}".format(self.distance)
		
		else:
			return "("+",".join([str(x) for x in self.children])+"):"+"{:.7f}".format(self.distance)
		
	def push(self, name = "", distance =0.0, index = 0):
		
		newchild = node(self,name,distance,index)
		
		return newchild
	
	def str(self):
		
		return str(self)
	
	def build(self,text):
		
		if text == "":
			
			return self
		
		elif text[-1] == ";":
			
			text = text[:-1]
			
		current_node = self
		allnames = []
		
		index = 0
		current_name = ""
		for char in text:
			
			if char in [" ", "\'"]:
				
				continue
			
			if char == "(":
				
				current_node = current_node.push()
				
			elif char == ")": 
				
				name,distance = (current_node.name.split(":")+["0.0"])[:2]
				
				name = name.split()[0].split(" ")[0].split("\t")[0]
				
				if len(name)>0:
					
					current_node.name = name
					
					allnames.append(name)
					
					if len(current_node.children) == 0:
						current_node.index = index
						index += 1
						
				else:
					
					allnames.append(current_node.name)
					
				try:
					current_node.distance = float(distance.strip())
					
				except:
					
					current_node.distance = 0.0
					
				current_node = current_node.parent
				
				if len(current_node.name) == 0 :
					
					current_node.name = "bh_"+str(len(allnames))
					
				else:
					print(current_node.name)
					
					
			elif char == "," or char ==";":
				
				name,distance = (current_node.name.split(":")+["0.0"])[:2]
				
				name = name.split()[0].split(" ")[0].split("\t")[0]
				
				if len(name)>0:
					
					current_node.name = name
					
					allnames.append(name)
					
					if len(current_node.children) == 0:
						current_node.index = index
						index += 1
						
				else:
					
					allnames.append(current_node.name)
					
				try:
					current_node.distance = float(distance.strip())
					
				except:
					
					current_node.distance = 0.0
					
					
				if current_node.parent is not None:
					current_node = current_node.parent.push()
					
			else:
				
				current_node.name += char
				
				
		return self
	
	
	def combineclades(self, dm , combinedindex, cutoff = 0.00):
		
		if len(self.children) == 0:
			return self
		
		
		self.children[0] = self.children[0].combineclades(dm , combinedindex, cutoff)
		self.children[1] = self.children[1].combineclades(dm , combinedindex, cutoff)
		
		if len(self.children[0].children) > 1 or len(self.children[1].children) > 1 :
			
			return self
		
		lindex = self.children[0].index
		rindex = self.children[1].index
		
		if self.children[0].distance > 0.0 or self.children[1].distance > 0.0 or dm[lindex, rindex] < max(dm[rindex, rindex],dm[lindex, lindex]) :
			
			return self
		
		self.children[0].distance = self.distance
		self.children[0].parent = self.parent
		
		combinedindex.add(rindex)
		
		
		return self.children[0]
	
	def filterindex(self,filterlist):
		
		if len(self.children) == 0:
			
			if self.index in filterlist:
				return None
			
			return self
		
		ldistance = self.children[0].distance
		rdistance = self.children[1].distance
		
		self.children[0] = self.children[0].filterindex(filterlist)
		self.children[1] = self.children[1].filterindex(filterlist)
		
		if self.children[0] is None and self.children[1] is None:
			
			return None
		
		elif self.children[1] is None:
			
			self.children[0].distance += self.distance
			self.children[0].parent = self.parent
			
			return self.children[0]
		
		elif self.children[0] is None:
			
			self.children[1].distance += self.distance
			self.children[1].parent = self.parent
			
			return self.children[1]
		
		return self



def overlapregion(segments):
	
	allcoordi = [x for y in segments for x in y]
	
	allcoordi_sortindex = sorted(list(range(len(allcoordi))), key = lambda x: allcoordi[x])
	
	merged = []
	cover = 0
	for i,x in enumerate(allcoordi_sortindex ):
		
		if i % 2:
			cover -= 1
			if cover == 0:
				merged.append(str(allcoordi[x]))
		else:
			cover += 1
			if cover == 1:
				merged.append(str(allcoordi[x]))
				
	return "_".join(merged)


def pathinformation(alignfile):
	
	allpathaligns = dict()
	with open(alignfile, mode = 'r') as f:
		
		for line in f:
			if len(line):
				line = line.strip().split()
				if len(line) < 3:
					continue

				name, pathes, alignregion_q, alignregion_r = line[0],line[1].replace("<",">").split(">")[1:], line[3].split(";"), line[4].split(";")
				
				pathaligns = cl.defaultdict(list)
				for path, region_q, region_r in zip(pathes, alignregion_q, alignregion_r):
					pathaligns[path].append(region_q+"_"+region_r)
					
				allpathaligns[name] = "&".join(["{}:{}".format(name,  "~".join(data)) for name, data in pathaligns.items()])

	return allpathaligns

def getgroupkmercounts(KmerReader, groups):
	
	groups_dict = dict()
	for i,group in enumerate(groups):
		for x in group:
			groups_dict[x] = i
			
			
	groups_kmercounts = [0 for x in groups]
	
	for kmerindex in range(len(KmerReader.kmerslist)):
		
		counts = list(KmerReader.matrix[kmerindex])[:KmerReader.rowlength[kmerindex]]
		
		counts_groups  = set([groups_dict[x] for x in counts])
		
		for groupindex in counts_groups:
			groups_kmercounts[groupindex] += 1
			
	return groups_kmercounts

	
def main(args):


	excludenames = set([])	
	if len(args.exclude):
		excludenames = set(args.exclude.split(","))
	
	allpathaligns = cl.defaultdict(str)
	if len(args.align) == 0 and os.path.isfile(args.seq + "allgraphalign.out"):
		args.align = args.seq + "allgraphalign.out"

	if len(args.align):
		allpathaligns = pathinformation(args.align)
	
	outfile = args.output
	
	name = "_".join(args.seq.split("/")[-1].split("_")[:2])
	
	kfilesize = os.stat(args.kmer).st_size
	
	KmerReader = KmerData(args.seq, msize = kfilesize)
	pathtitles = KmerReader.LoadKmers(args.kmer)
	
	genenum = len(KmerReader.sampleslist)
	kmernum = len(KmerReader.kmerslist)
	
	useindex = set([i for i,name in enumerate(KmerReader.headers) if name.split()[1].split(":")[0].split("#")[0] in excludenames])

	if os.path.isfile(args.seq +"_tree.ph") and os.stat(args.seq +"_tree.ph").st_size > 1 and ( ( os.path.isfile(args.seq +"_norm.txt") and os.stat(args.seq +"_norm.txt").st_size > 1) or ( os.path.isfile(args.seq +"_norm.gz") and os.stat(args.seq +"_norm.gz").st_size > 1) ):
		
		with open(args.seq +"_tree.ph", mode = 'r') as f:
			treetext = f.read().strip()
		
		if  ( os.path.isfile(args.seq +"_norm.txt") and os.stat(args.seq +"_norm.txt").st_size > 1): 
			with open(args.seq+ "_norm.txt", mode = 'r') as f:
				sqmatrix = []
				for line in f:
					if len(line) == 0:
						continue
				
					sqmatrix += list(map(float,line.strip().replace(","," ").split(" ")))
		else:	
			with gzip.open(args.seq+ "_norm.gz", mode = 'rt', encoding='utf-8') as f:
				sqmatrix = []
				for line in f:
					if len(line) == 0:
						continue
					sqmatrix += list(map(float,line.replace(","," ").strip().split(" ")))
		sqlength = len(sqmatrix)
		if sqlength != (int(sqlength**0.5))**2:
			sqmatrix = MatrixUppertoFull(sqmatrix, int((2*sqlength)**0.5) )

	else:
		KmerReader.initiatematrix()
		KmerReader.ReadKmers(args.seq)
	
		sqmatrix = KmerReader.matrix.SquareMatrix(0)
		
		treetext = projectiontree(np.reshape(sqmatrix, (genenum, genenum)), KmerReader.sampleslist)
		
	filterindex = set([])
	
	if args.ifmerge:

		treetext = node().build(treetext).combineclades(np.reshape(sqmatrix, (genenum, genenum)), filterindex).str() + ";"

		filterindex.update(useindex)
	
	if len(filterindex):

		treetext = node().build(treetext).filterindex(useindex).str() + ";"

		sqmatrix = MatrixFilter(sqmatrix, genenum, filterindex)

		genenum -= len(filterindex)

	#np.savetxt(args.output + "_norm.txt", np.reshape(sqmatrix, (len(KmerReader.sampleslist), len(KmerReader.sampleslist)) ), delimiter=',', fmt='%1.4f')


	groups = ""

	"""
	if os.path.isfile(args.seq +"_groups.txt"):
		with open(args.seq +"_groups.txt", mode = 'r') as f:
			groups = f.read().splitlines()
	"""

	groups,sample_groupindex = group_matrix(sqmatrix, genenum, 0.8)
	groupkmercounts = getgroupkmercounts(KmerReader,groups)
	groups = ["@"+str(int(c))+"\t"+",".join([str(x) for x in sorted(group)])+"," for c, group in zip(groupkmercounts,groups) if len(group) ]	
	#groups = [",".join([str(x) for x in group])+"," for group in groups if len(group) ] 

	groups = "\n".join([ "@"+x for x in groups])
	
	sqmatrix = MatrixFulltoUpper(sqmatrix, genenum)
	
	sqmatrix = "\n$".join([  " ".join(list(map(str,sqmatrix[i]))) for i in range(genenum)])
	
	kfilesize = os.stat(args.kmer).st_size
	
	filternames = [x for i,x in enumerate(KmerReader.sampleslist) if i in filterindex]
	keepnames = [x for i,x in enumerate(KmerReader.sampleslist) if i not in filterindex]

	KmerReader.initiatematrix()	
	KmerReader.ReadKmers(args.seq, filternames)
	
	genekmercounts = " ".join([str(x) for i,x in enumerate(KmerReader.genekmercounts) ])

	titletext = [""]
	for i,x in enumerate(KmerReader.headers):
		x = x.strip().split("\t")[:2]
		x = "\t".join(x +[allpathaligns.get(x[0][1:],'NA:NA')])

		if i in useindex:
			continue
		elif i in filterindex:
			titletext[-1]+=";"+x
		else:
			titletext.append(x)
	titletext = "\n".join(titletext[1:])

	
	with open(outfile, mode = 'w') as f:
		
		f.write("#"+"\t".join([name, str(kmernum), str(genenum), " "*30])+"\n")
		f.write("%"+genekmercounts + "\n")
		f.write("!"+treetext +"\n")
		if len(pathtitles):
			f.write(pathtitles + "\n")
		f.write(titletext+"\n")

		#f.write("".join(["\n>"+x[1:] if i not in filterindex else ";>"+x[1:] for i,x in enumerate(KmerReader.headers) if i not in exludeindex ])+"\n" )
		#f.write("L"+" ".join([str(KmerReader.pathsize[i]) for i in range(KmerReader.numpath + 1)])+"\n")
		#f.write("$"+sqmatrix+"\n")
		f.write(groups+"\n")
		if len(KmerReader.kmerslist) == 0:
			return
		
		allkmers = sorted( list(KmerReader.kmerslist.keys()) , key = lambda x: KmerReader.kmerslist[x])
		
		allkmers_sortindex = sorted(range(len(allkmers)), key = lambda x: hashrow(x, genenum, KmerReader ))
		
		#allkmers_sort = sorted(list(KmerReader.kmerslist.keys()), key = lambda x: KmerReader.kmerslist[x])
		
		lastcounts = None
		lastsign = ''
		for kmerindex in allkmers_sortindex:
			
			kmer = allkmers[kmerindex]
			
			counts = list(KmerReader.matrix.getkmerrow(kmerindex))
			
			counts = sorted(counts)
			
			kmerlistindex, kmerpath, kmerloc, kmerflag, kmerratio = KmerReader.kmerslist[kmer]
	
			theset = set(counts)
			
			if len(theset) <= len(keepnames)*0.5 or len(counts) != len(theset):
				sign = "+"
				row = ",".join([intencode(x) for x in counts]+[""])
			else:
				sign = "-"
				row = ",".join([intencode(x) for x in range(len(keepnames)) if x not in theset]+[""])
			

	
			if sign == lastsign and counts == lastcounts:
				line = ["&" + ( "_" if sign == '-' else "=" ) , intencode(max(0,int(kmerflag))) + "|" + intencode( min(1000, int(100 * float(kmerratio) + 0.5)) ) , intencode(max(0,int(kmerpath)), length = 3) , intencode(max(0, int(kmerloc)), length = 5), intencode(kmer, length = 11), ""]
			else:   
				line = ["&" + sign, intencode(max(0,int(kmerflag))) + "|" + intencode( min(1000, int(100 * float(kmerratio) + 0.5)) )  ,intencode(max(0,int(kmerpath)), length = 3) , intencode(max(0, int(kmerloc)), length = 5), intencode(kmer, length = 11), row]
				
			f.write("\t".join(line)+"\n")
			
			lastcounts = counts
			lastsign = sign
			
			
def run():
	"""
		Parse arguments and run
	"""
	parser = argparse.ArgumentParser(description="program compiles sequences to kmer matrix")
	parser.add_argument("-s", "--seq", help="path to input data file", dest="seq", type=str, required = True)
	parser.add_argument("-k", "--kmer", help="path to input data file", dest="kmer", type=str, required = True)
	parser.add_argument("-o", "--output", help="path to output file", dest="output", type=str, required = True)
	parser.add_argument("-a", "--align", help="path to output file", dest="align", type=str, default = "")
	parser.add_argument("-m", "--ifmerge", help="path to output file", dest="ifmerge", type=int, default = 1)
	parser.add_argument("-e", "--exclude", help="path to output file", dest="exclude", type=str, default = "")
	parser.set_defaults(func=main)
	args = parser.parse_args()
	args.func(args)
	
	
if __name__ == "__main__":
	run()
