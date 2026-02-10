#!/usr/bin/env python3
import argparse
import sys
import collections as cl
def main(args):
    
    ifdup = 0
    iffilter = 0
    headers = set()
   
    prefix = args.input.split("/")[-1].split("_")[0]
 
    w = open(args.output, mode = 'w')
    with open(args.input, mode = 'r') as f:
        
        for line in f:
            if len(line) ==0:
                continue
            if line[0] == ">":
                name = "_".join(line.split()[0].split("_")[1:])
                iffilter = 0
                if name in headers:
                    iffilter =1 
                    ifdup = 1
                headers.add(name)

                line = line.strip().split("_")
                line = "_".join([">"+prefix]+line[1:])+"\t"+line[0][1:].replace("N","_") + "\n"
            if not iffilter:
                w.write(line)
        
    w.close()
    
    if ifdup:
        print(args.input)
                


def run():
    """
        Parse arguments and run
    """
    parser = argparse.ArgumentParser(description="program distract overlap genes and alignments on contigs")
    parser.add_argument("-i", "--input", help="path to input data file",dest="input", type=str, required=True)
    parser.add_argument("-o", "--output", help="path to input data file",dest="output", type=str, required=True)
    parser.set_defaults(func=main)
    args = parser.parse_args()
    args.func(args)
    
    
if __name__ == "__main__":
    run()
