#!/usr/bin/python

import collections as cl
import math
import numpy as np
import argparse
import os
import gzip
script_folder = os.path.dirname(os.path.abspath(__file__))

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
                                
                        firstrow = [float(x) for x in firstrow]
                        
                        self.dim = len(firstrow)
                        
                        self.data = np.zeros(self.dim*(self.dim+1)//2, dtype=np.int32)
                        
                        self.data[0:self.dim] = firstrow 
                        
                        index = self.dim
                        for row in f:
                                
                                row = row.decode().strip().strip(",").split(",")[:-1]
                                self.data[index:(index+len(row))] = [float(x) for x in row]
                                index += len(row)
                                        
                return self

def readfasta(reffile):

        with open(reffile, mode = 'r') as f:

                readtext = f.read()

                readtext = readtext.split(">")[1:]

                #readtext = [read for read in readtext if samplename not in read.splitlines()[0]]

                headers = [read.splitlines()[0] for i,read in enumerate(readtext)]
                text = ["".join(read.splitlines()[1:]) for i,read in enumerate(readtext)]
                names = [header.split()[1].split("#")[0] for header in headers]

        return text, headers, names 

def group_matrix(matrix, size, similar):

        index = 0
        vars = []
        for i in range(size):
                vars.append(matrix[index])
                index += (size - i)
        
        groups = [[] for x in range(size)]
        group_edges = [x for x in range(size)]
        for i in range(size):
                
                indexstart = i*size - (i*(i+1)//2)
                var1 = vars[i] 
                for j in range(i+1, size):
                        
                        edge =float( matrix[indexstart + j ])
                        if edge* edge > similar * similar * float(var1) * float(vars[j])  :
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



def readfasta(reffile):

        with open(reffile, mode = 'r') as f:

                readtext = f.read()

                readtext = readtext.split(">")[1:]

                #readtext = [read for read in readtext if samplename not in read.splitlines()[0]]

                headers = [read.splitlines()[0] for i,read in enumerate(readtext)]
                text = ["".join(read.splitlines()[1:]) for i,read in enumerate(readtext)]
                names = [header.split()[1].split("#")[0] for header in headers]

        return text, headers, names 


def makereverse(seq):

        tran=str.maketrans('ATCGatcg', 'TAGCtagc')

        return seq[::-1].translate(tran)



def similar(table1, table2):


        if len(table1)>len(table2):
                larger, smaller = table1, table2
        else:
                larger, smaller = table2, table1

        total_counts = sum(list(smaller.values()))


        concordance = 0
        for kmer, count in smaller.items():

                larger_count = larger[kmer]

                concordance += min(count, larger_count)

        return concordance/total_counts


class KmerData:

        def __init__(self,samplelist, selectlist, matrix):

                self.selectkmers = selectlist
                self.samplelist = samplelist

                self.kmerlist = dict()
                self.tran=str.maketrans('ATCGatcg', 'TAGCtagc') 
                self.matrix = matrix
                self.matchscores = 0.0  
                self.weightscores = 0.0

        def getkmerindex(self, kmer):

                kmer = max(kmer, kmer[::-1].translate(self.tran))

                return self.kmerlist.get(kmer, -1)       

        def getsampleindex(self,kmer):

                kmerindex = self.getkmerindex(kmer)

                if kmerindex != -1:

                        samples = set(self.matrix.getkmersamples(kmerindex).keys())

                else:
                        samples = set()

                return samples 

        def countkmer(self, seq, index ,kmertype = 31):

                kmer_counter = 0
                for posi in range(len(seq)-kmertype+1):

                        kmer = seq[posi:posi+kmertype]  

                        kmer = max(kmer.upper(), kmer[::-1].translate(self.tran).upper())

                        if len(self.selectkmers) and kmer not in self.selectkmers:
                                continue

                        kmerindex = self.kmerlist.get(kmer, -1) 

                        if kmerindex < 0:

                                kmerindex = len(self.kmerlist)
                                self.kmerlist[kmer] = kmerindex
                                self.matrix.addkmer()

                        self.matrix.add(kmerindex, index)
                        kmer_counter += 1

                return kmer_counter


def main(args):


        inputfile = args.input

        outputfile = args.output

        kmerfile = args.kmer

        selectlist = set()

        #os.system(f"{script_folder}/kmernorm -i {inputfile} -k {kmerfile} -o {output}_norm.gz")
        
        matrix = Matrix().load(outputfile+"_norm.gz")

        seqs, headers, names = readfasta(inputfile)

        groups, sample_groupindex = group_matrix(matrix.data, matrix.dim, args.similar)

        for group_index, group in enumerate(groups):

                if len(group) == 0:
                        continue

                with open(outputfile+"{}.fa".format(group_index+1), mode = "w") as f:

                        for index in group:

                                f.write(">{}\n{}\n".format(headers[index],seqs[index]))



def run():
        """
                Parse arguments and run
        """
        parser = argparse.ArgumentParser(description="program determine psuedogene")
        parser.add_argument("-i", "--input", help="path to input data file", dest="input", type=str, required=True)
        parser.add_argument("-o", "--output", help="path to output file", dest="output", type=str, required=True)
        parser.add_argument("-k", "--kmer", help="path to kmer file", dest="kmer", type=str, default="")
        parser.add_argument("-s", "--similar", help="path to kmer file", dest="similar",type=float, default=0.5)
                        
        parser.set_defaults(func=main)
        args = parser.parse_args()
        args.func(args)


if __name__ == "__main__":
        run()
