#!/usr/bin/env python3


import collections as cl
import math
import numpy as np
import scipy 
import pandas as pd
import argparse

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



class Kmer_annotation:
    
        def __init__(self,  kmersize =31):
            
                self.tran=str.maketrans('ATCGatcg', 'TAGCtagc') 
                self.kmerlist = set([])
                self.kmersign = cl.defaultdict(lambda x: 1)
                self.kmersize = kmersize
                self.intoperator = 4**(kmersize-1)
                self.contigsign = cl.defaultdict(int)
            
    
        def annotate(self, seq, totalsign = 0):
            
                kmersize = self.kmersize
                foundposi = []
            
                self.current_k = 0
                self.reverse_k = 0
                self.current_size = 0
            
                lastposi = -10000000
            
            
                for posi,char in enumerate(seq):
                    
                        char = char.upper()
                    
                        if char not in ['A','T','C','G']:
                            
                                self.current_k = 0
                                self.reverse_k = 0
                                self.current_size = 0
                                continue
                    
                        if self.current_size >= kmersize:
                            
                                self.current_k %= self.intoperator
                                self.current_k <<= 2
                                self.current_k += ['A','C','G','T'].index(char)
                            
                                self.reverse_k >>= 2
                                self.reverse_k += (3-['A','C','G','T'].index(char)) * self.intoperator
                            
                        else:
                            
                            
                                self.current_k <<= 2
                                self.current_k += ['A','C','G','T'].index(char)
                                self.reverse_k += (3-['A','C','G','T'].index(char)) * ( 1 << (2*self.current_size))
                            
                                self.current_size += 1
                            
                        if self.current_size >= kmersize:  
                            
                                kmer = max(self.current_k, self.reverse_k)
                            
                                if kmer == self.current_k:
                                        thesign1 = 1
                                else:
                                        thesign1 = -1
                                    
                                thesign2 = self.kmersign.get(kmer, 0)
                                if thesign2 == 0:
                                    
                                        pass
                                    
                                else:
                                    
                                        totalsign += thesign2 * thesign1
                                    
                if totalsign >= 0:
                        totalsign = max(1, totalsign)
                
                    
                for posi,char in enumerate(seq):
                    
                        char = char.upper()
                    
                        if char not in ['A','T','C','G']:
                            
                                self.current_k = 0
                                self.reverse_k = 0
                                self.current_size = 0
                                continue
                    
                        if self.current_size >= kmersize:
                            
                                self.current_k %= self.intoperator
                                self.current_k <<= 2
                                self.current_k += ['A','C','G','T'].index(char)
                            
                                self.reverse_k >>= 2
                                self.reverse_k += (3-['A','C','G','T'].index(char)) * self.intoperator
                            
                        else:
                            
                            
                                self.current_k <<= 2
                                self.current_k += ['A','C','G','T'].index(char)
                                self.reverse_k += (3-['A','C','G','T'].index(char)) * ( 1 << (2*self.current_size))
                            
                                self.current_size += 1
                            
                        if self.current_size >= kmersize:  
                            
                                kmer = max(self.current_k, self.reverse_k)
                            
                                if kmer == self.current_k:
                                        thesign1 = 1
                                else:
                                        thesign1 = -1
                                    
                                thesign2 = self.kmersign.get(kmer, 0)
                                if thesign2 == 0:
                                    
                                        self.kmersign[kmer] = thesign1 * np.sign(totalsign)
                                    
                                        totalsign += 1
                                    
                                
                                    
                                    
                                    
                return totalsign
    
    
    
def main(args):
    
    
        infile = args.input
    
        outfile = args.output
    
    
        annotator = Kmer_annotation()
    
    
        with open(infile, mode = 'r') as f:
            
                reads = f.read().split(">")[1:]
                reads = [read.splitlines() for read in reads]
                names = [read[0].split()[0] for read in reads]
                reads = {read[0].split()[0]: (read[0],"".join(read[1:])) for read in reads}
          
        allnames = names 
        contigsize = cl.defaultdict(int)
    
        nametocontig = dict()
        for name, (header,seq) in reads.items():
            
            contig = header.split()[1].split(":")[0] 
            size = len(seq)
            contigsize[contig] += size
            nametocontig[name] = contig
            
        names_sort = sorted(names, key = lambda x: ( contigsize[nametocontig[x]], len(reads[x][1]) ), reverse = 1)
        nametosign = dict()
        contigsign = cl.defaultdict(int)
        for name in names_sort:
            
                (header, seq) = reads[name]
                
                sign = annotator.annotate(seq, min(500, contigsign[nametocontig[name]]//100))
                contigsign[nametocontig[name]] += sign 
            
                nametosign[name] = sign 
            
        with open(outfile, mode = 'w') as w:
            
                for name in allnames:
                    
                        (header, seq) = reads[name]
                        thesign = nametosign[name]
                        header = header.split()
                        oldstrd = ""
                        if header[1][-1] in ['+','-']:
                                pass
                                oldstrd = header[1][-1]
                                header[1] = header[1][:-1]
                                
                        if thesign  < 0:
                                seq = makereverse(seq)
                                strd = "-" if oldstrd != '-' else '+'
                        else:
                                strd = "+" if oldstrd != '-' else '-' 

                        header[1] += strd                          
 
                        header = "\t".join(header)
                    
                        w.write(">"+header + "\n"+ seq+"\n")
                    
                    
def run():
        """
                Parse arguments and run
        """
        parser = argparse.ArgumentParser(description="program determine strand")
        parser.add_argument("-i", "--input", help="path to input data file", dest="input", type=str, required=True)
        parser.add_argument("-o", "--output", help="path to output file", dest="output", type=str, required=True)
        parser.set_defaults(func=main)
        args = parser.parse_args()
        args.func(args)
    
    
if __name__ == "__main__":
        run()
