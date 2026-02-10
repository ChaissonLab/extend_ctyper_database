#!/usr/bin/env python3


import os 
import argparse
import sys
import numpy as np
import math
import collections as cl
import gzip

mainchr = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM']

def normtotdm(matrix):
    
    
    matrix = np.matrix(matrix)
    
    dm = np.ones((matrix.shape))
    vars = [matrix[i,i] for i in range(len(matrix))]
    for i in range(len(matrix)):
        
        var1 = max(1.0,vars[i])
        for j in range(i+1):
            
            var2 = max(1.0,vars[j])
            
            dm [i,j] = 1 - matrix[i,j]/math.sqrt((var1*var2))
            dm [j,i] =  dm [i,j]
            
    return dm

def load_upper_tri_norm(path: str) -> np.ndarray:
    # reads an upper-triangular CSV (often ragged) and returns full symmetric matrix
    opener = gzip.open if path.endswith(".gz") else open
    rows = []
    with opener(path, "rt") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split(",")
            if parts and parts[-1] == "":  # drop trailing comma
                parts = parts[:-1]
            # keep non-empty tokens only (in case there are blanks)
            vals = [float(x) for x in parts if x != ""]
            rows.append(vals)
            
    dim = len(rows[0])
    M = np.zeros((dim, dim), dtype=np.float32)
    
    for i, vals in enumerate(rows):
        # Expect upper-tri ragged: len(vals) == dim - i
        if len(vals) == dim - i:
            M[i, i:dim] = vals
        else:
            # fallback: if it's already full width, just take row i
            M[i, :len(vals)] = vals
            
    # mirror to lower triangle
    M = M + M.T - np.diag(np.diag(M))
    return M

def main(args):
    
    
    inputfile = args.input
    outputfile = args.output
    
    normfile = args.input + "_norm.txt.gz"
    treefile = args.input + "_tree.ph"
    
    normmatrix = load_upper_tri_norm(normfile)
    
    dm = normtotdm(normmatrix)
    
    header_info = dict()
    headers = dict()
    with open(args.input, mode = 'r') as f:
        for line in f:
            if len(line) ==0 or line[0] != ">":
                continue
            lines = line.strip().split()
            headers[lines[0][1:]] = lines[1]
            header_info[lines[0][1:]] = lines[-1]
    nonrefindexes = []
    refindexes = []
    for i, header in enumerate(headers.keys()):
        
        if len(args.ref):
            if args.ref in header:
                refindexes.append(i)
                nonrefindexes.append(i)
            else:
                nonrefindexes.append(i)
        else:
            loci = headers[header]
            if loci.split(":")[0].startswith(mainchr):
                refindexes.append(i)
                nonrefindexes.append(i)
            else:
                nonrefindexes.append(i)
                
    headers_list = list(headers.keys())
    
    allpairs = []
    for i in nonrefindexes:
        
        for j in refindexes:
            
            allpairs.append( ( dm[i,j], i,j) )
            
    allpairs = sorted(allpairs)
    
    matched_query = set()
    matched_pair = set()
    best_matches = cl.defaultdict(list)
    all_matches = cl.defaultdict(list)
    
    for (distance, i, j) in allpairs:
        
        haplo1 = "_".join(headers_list[i].split('_')[:-1])
        haplo2 = headers_list[j]
        
        thepair = (min(haplo1,haplo2), max(haplo1,haplo2))
        
        if distance > 0.5:
            break
                
        ifdup = 0
        if thepair in matched_pair:
            ifdup = 1
        
        elif headers_list[i] in matched_query:
            ifdup = 2
            
        best_matches[headers_list[i]] = [ (distance, haplo2 , ifdup  ) ]
        if ifdup != 1:
            matched_pair.add(thepair)
            matched_query.add(headers_list[i])
            all_matches[headers_list[i]].append( (distance, headers_list[j], ifdup) )
            
        elif distance < 0.1:
            all_matches[headers_list[i]].append( (distance, headers_list[j], ifdup) )
            
    tags = ["Pri","Dup", "Cov"]
    
    with open(outputfile, mode = 'w') as f:
        
        for header in headers:
            if header in best_matches:
                
                for (distance, ref, ifdup) in best_matches[header][:1]:
                    matches =  ";".join([f"{x[1]}:{round(x[0],6)}" for x in sorted(all_matches[header])])
                    
                    f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(header, ref, headers[header], headers[ref], int("ENSE" in header_info[header]) , distance,tags[ifdup],matches))
                    
                    break
                
            else:
                if args.ref not in header: 
                    f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(header, "NA", headers[header], "NA", int("ENSE" in header_info[header]),0.0,"Novel",""))
                else:
                    f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(header, "NA", headers[header], "NA", int("ENSE" in header_info[header]),0.0,"Ref",""))
                    
                    
                    
    return 0

def run():
    """
        Parse arguments and run
    """
    parser = argparse.ArgumentParser(description="program determine genes")
    parser.add_argument("-i", "--input", help="path to input data file", dest="input", type=str, required=True)
    parser.add_argument("-o", "--output", help="path to output file", dest="output", type=str, required=True)
    parser.add_argument("-r", "--ref", help="path to output file", dest="ref", type=str, default="CHM13_h1")
    parser.set_defaults(func=main)
    args = parser.parse_args()
    args.func(args)
    
    
if __name__ == "__main__":
    run()
