#!/usr/bin/env python3


import os
import argparse
import collections as cl
import pandas as pd
from bitarray import bitarray

def makereverse(sequence):
    
    tran=str.maketrans('ATCGatcg', 'TAGCtagc')
    
    return sequence[::-1].translate(tran)

def kmerencode(kmer):
    
    kmerint = 0
    for base in kmer:
        
        kmerint *= 4
        if base=='a' or base =='A':
            
            kmerint += 0
            
        elif base=='t' or base =='T':
            
            kmerint += 3
            
        elif base=='c' or base =='C':
            
            kmerint += 1
        elif base=='g' or base =='G':
            
            kmerint += 2
            
    return kmerint

def kmerdecode(kmerint, size = 31):
    
    index = 0
    text = ['A' for a in range(size)]
    
    while kmerint:
        
        text[index] = "ACGT"[kmerint % 4]
        kmerint //= 4
        index += 1
        
    return "".join(text[::-1])


                    
class KmerReader:
    
    def __init__(self,  kmersize = 31):
        
        self.kmersize = kmersize
        self.kmerindexes = cl.defaultdict(int)
        
        
    def loadsamples(self, infile):
        
        self.samplesizes = [0]
        self.samples = []
        with open(infile, mode = 'r') as f:
            
            for line in f:
                
                if line.startswith(">"):
                    
                    name = line.split()[0][1:]
                    self.samples.append(name)
                    self.samplesizes.append(0)
                else:
                    self.samplesizes[-1] += len(line.strip())
                    
        self.samplesizes = self.samplesizes[1:]
        self.samplenum = len(self.samples)
        return self
    
    def loadkmers(self, kmerfile):
        
        kmercount = 0
        with open(kmerfile, mode = 'r') as f:
            
            for line in f:
                
                if line.startswith(">") == False:
                    
                    kmercount += 1
                    self.kmerindexes[kmerencode(line.strip().split()[0])] = kmercount
                    
        self.kmernum = kmercount + 1
        
        
        return self
    
    def read(self, file):
        
        kmercounter = [0] * (self.kmernum + 1)
        
        intoperator = 4**(self.kmersize-1)
        kmersize = self.kmersize
        
        current_k = 0
        reverse_k = 0
        current_size = 0
        groupindex =-1 
        
        currentline = ""
        with open(file,mode = 'r') as r:
            
            for line in r:
                
                line = line.strip()
                
                if line.startswith(">"):
                    
                    current_k=0
                    reverse_k=0
                    current_size = 0
                    
                    continue
                
                for char in line.upper():
                    
                    if char not in ['A','C','G','T']:
                        
                        current_k = 0
                        reverse_k = 0
                        current_size = 0
                        continue 
                    
                    kmercode = ['A','C','G','T'].index(char)
                    
                    if current_size >= self.kmersize:
                        
                        current_k %= intoperator
                        current_k <<= 2
                        current_k += kmercode
                        
                        reverse_k >>= 2
                        reverse_k += (3-kmercode) * intoperator
                        
                    else:
                        
                        
                        current_k <<= 2
                        current_k += kmercode
                        reverse_k += (3-kmercode) * ( 1 << (2*current_size))
                        
                        current_size += 1
                        
                    kmer = max(current_k, reverse_k)
                    
                    if current_size == kmersize:
                        
                        kindex = self.kmerindexes.get(kmer, 0)
                        
                        if kindex > 0:
                            kmercounter[kindex] += 1
        return kmercounter
    
    
    
def selectkmers(reader, count1 ,count2, outputfile):
    
    kmertable = reader.kmerindexes

    with open(outputfile, mode = 'w') as w:
        
        for kmer,kindex in reader.kmerindexes.items():
                        
            if count1[kindex] == count2[kindex] :
                
                w.write(">\n{}\n".format(kmerdecode(kmer)))
                
                
                
def main(args):
    
    kmertype = 31
    
    kmerfile = args.kmer
    
    infile1 = args.f1
    infile2 = args.f2
    
    kmeroutfile = args.output
    
    
    reader = KmerReader().loadkmers(kmerfile)
    count1 = reader.read(infile1)
    count2 = reader.read(infile2)
        
    selectkmers(reader, count1 ,count2, kmeroutfile)
    
    
    
def run():
    """
        Parse arguments and run
    """
    parser = argparse.ArgumentParser(description="program determine psuedogene")
    parser.add_argument("-f1", "--file1", help="path to input data file",dest="f1", type=str, required=True)
    parser.add_argument("-f2", "--file2", help="path to input data file",dest="f2", type=str, required=True)
    parser.add_argument("-k", "--kmer", help="path to input data file",dest="kmer", type=str, required=True)
    parser.add_argument("-o", "--output", help="path to input data file",dest="output", type=str, required=True)
    parser.set_defaults(func=main)
    args = parser.parse_args()
    args.func(args)
    
    
if __name__ == "__main__":
    run()
