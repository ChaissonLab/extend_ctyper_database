#!/usr/bin/env python3

import multiprocessing as mul
import argparse
import os
import subprocess
import datetime
from graphfilesplit import filesplit
from graphinserts import insertdistract
from graphcigar import graphcigar
from graphcleanrepeats import cleanrepeats
from graphtolinear import graphtolinear
import re

script_folder = os.path.dirname(os.path.abspath(__file__))
lock1 = mul.Lock()

ifblast_g = 1


def checkmask(fasta_file):
    
    bash_cmd = f'''
    grep -v '^>' "{fasta_file}" | tr -d '\\n' | tr -cd 'a-z' | wc -c
    '''
    result = subprocess.run(["bash", "-c", bash_cmd], capture_output=True, text=True)
    lowercase_count = int(result.stdout.strip())
    
    return lowercase_count

def fileordered(allfiles, priors=None):
   
    def file_has_Nn(fa_path):
            with open(fa_path, "r") as f:
                for line in f:
                    if line.startswith('>'):   # header
                        continue
                    if ('N' in line) or ('n' in line):
                        return True
            return False
    
    def header_has_chr_colon(fa_path):
        _chr_pat = re.compile(r"^chr(?:[0-9]+|X|Y):")
        with open(fa_path, "r") as f:
            return _chr_pat.match(f.readline().strip()) is not None
        
    priors = priors or set()
    
    def score(path):
        fn = os.path.basename(path)
        stem = fn[:-3] if fn.endswith(".fa") else os.path.splitext(fn)[0]
        sz = os.path.getsize(path)
        
        if stem in priors:
            return 10e20 + sz
        elif ("NC_0609" in fn) or ("CHM13_" in fn):
            return 5e20 + sz
        elif "HG002_" in fn:
            return 10e17 + sz
        if file_has_Nn(path):
            sz = sz - 10e20
        if "HG38alt" in fn:
            return 10e18 + sz
        elif (("_chr" in fn) or ("HG38_" in fn)):
            return (10e19 + sz) if header_has_chr_colon(path) else (10e18 + sz)
        
        else:
            return sz
        
    return sorted(allfiles, key=score, reverse=True)

def selfalign(input,output,nthreads, ifblast = 0):
    
    if ifblast:
        dbcmd = "bash {}/runmakeblastdb -in {} -out {}_db ".format(script_folder,input, input)
        os.system(dbcmd)
        cmd = "bash {}/runblastn  -task megablast  -query {} -db {}_db -gapopen 10 -gapextend 2 -word_size 30  -perc_identity 95 -dust yes -lcase_masking -evalue 1e-200 -outfmt 17 -out {}_selfblast.out  -num_threads {} -max_target_seqs 100 ".format( script_folder,input , input, input,nthreads)
        os.system(cmd)
        totalclean = cleanrepeats(input, input+"_", nthreads)
        
        cmd = "{}/runwinnowmap.sh {} {} {} 0.95 150  > {}_selfblast.out".format(script_folder,  input+"_",input+"_", nthreads, input+"_")
        os.system(cmd)
        totalclean += cleanrepeats(input+"_", output, nthreads)
    else:  
        cmd = "{}/runwinnowmap.sh {} {} {} 0.95 150 > {}_selfblast.out".format(script_folder,  input,input, nthreads, input)
        os.system(cmd)
        
        totalclean = cleanrepeats(input, output, nthreads)

    return totalclean

def callblastn(query, graphfile,  output,nthreads, simi = 0.90, ifrepeat = 1):
    
    if ifrepeat:
        repeatopts = "-dust yes -lcase_masking"
    else:
        repeatopts = ""
        
    cmd = "bash {}/runblastn -task megablast -query {} -db {} -gapopen 10 -gapextend 2 -word_size 30   -perc_identity {} {} -evalue 1e-200  -outfmt 17 -out {}  -num_threads {} -max_target_seqs 100 ".format( script_folder, query,graphfile+"_db" , int(simi * 100),repeatopts,output,nthreads)
    os.system(cmd)
    
    
    
def callalign(query, graphfile,  output,nthreads, simi = 0.90,  ifhighqual = 1, ifrepeat = 1):
    if ifhighqual:
        callblastn(query, graphfile,  output,nthreads, simi, ifrepeat)
        
        
        #if ifrepeat:
            #repeatopts = "-dust yes -lcase_masking"
        #else:
            #repeatopts = ""
        
        #cmd = "bash {}/runblastn -task megablast -query {} -db {} -gapopen 10 -gapextend 2 -word_size 30   -perc_identity {} {} -evalue 1e-200  -outfmt 17 -out {}  -num_threads {} -max_target_seqs 100 ".format( script_folder, query,graphfile+"_db" , int(simi * 100),repeatopts,output,nthreads)
        #os.system(cmd)
        
        if (not os.path.exists(output) or os.path.getsize(output) <= 10):
            
            print(f"WARNING: alignment fail :{output}, realigning\n")
            
            callblastn(query, graphfile,  output,nthreads, simi, ifrepeat)
            
            #fallback_cmd = "bash {}/runblastn -task megablast -query {} -db {} -gapopen 10 -gapextend 2 -word_size 30 -perc_identity {} {} -evalue 1e-200 -outfmt 17 -out {} -num_threads {} -max_target_seqs 100".format(script_folder, query, graphfile+"_db", int(simi * 100), repeatopts, output, nthreads)
            
            #os.system(fallback_cmd)
            
            
        if ifrepeat:
            
            cmd = "{}/runwinnowmap.sh {} {} {} {} 150 > {}".format(script_folder,  query,graphfile, nthreads, simi, output+"_wm")
            os.system(cmd)
            if (not os.path.exists(output+"_wm") or os.path.getsize(output+"_wm") <= 10):
                os.system(cmd)
            os.system(f" ( cat  {output}_wm  >> {output} &&  rm {output}_wm  ) || true ")
            
    else:
        #cmd = "{}/runfast.sh {} {} {} 0.95 300 -c > {}".format(script_folder,  query,dbpath, nthreads,output)
        cmd = "{}/runwinnowmap.sh {} {} {} {} 150 > {}".format(script_folder,  query,graphfile, nthreads, simi, output)
        os.system(cmd)
        
        
def addnewseq(folder, graphfile, addfile,nthreads, simi, ifblast = 1):
    
    filealign = addfile + "_blastn.out"
    fileinserts = addfile + "_inserts.fasta"
    
    try:
        os.remove(fileinserts)
    except:
        pass
        
    if ifblast:
        callblastn( addfile,  folder+graphfile.split("/")[-1], filealign,nthreads, simi, 1)
    else:
        callalign( addfile,  folder+graphfile.split("/")[-1], filealign,nthreads, simi, 0)
        
    insertfiles = insertdistract( filealign, addfile, fileinserts, 1, 150)
    
    if len(insertfiles) and os.path.isfile(fileinserts) != False and os.path.getsize(fileinserts) > 10:
        os.system("cat {} >> {}".format(fileinserts, graphfile))
        
    else:
        return ""
    
    return insertfiles

def runcleanrepeats(folder, graphfile, newgraphfile,nthreads, ifhighqual = 1):


    size = 100000000
    
    ncycle = 1
    if ifhighqual:
        ncycle = 10
   
    theinput = graphfile
    theoutput = newgraphfile
    thetemp  = newgraphfile + ".temp"

    icycle = 0
    while icycle < ncycle and size > 300:
        if icycle > 0:
            ifhighqual = 0
            theinput = theoutput
            theoutput = newgraphfile if theoutput  == thetemp else thetemp
        size = selfalign(theinput, theoutput, nthreads, ifhighqual)
        icycle += 1

    if theoutput == thetemp:
        os.replace(thetemp, newgraphfile)
    elif icycle > 1:
        os.remove(thetemp)


    return newgraphfile

def editname(thefile, graphfile):
    
    size = 0
    with open(thefile, mode = 'r') as r:
        for line in r:
            if len(line) and line[0] != '>':
                size += len(line.strip())
                
                
    with open(thefile, mode = 'r') as r:
        with open(graphfile, mode = 'w') as w:
            for line in r:
                if len(line) and line[0] == '>':
                    
                    header = line.strip().split()
                    header[0] += "_{}_{}".format(0, size)
                    
                    if len(header[0]) > 50 -2:
                        header[0] = ">"+"_".join(header[0].split("_")[2:])
                    w.write("\t".join(header)+"\n")
                else:
                    w.write(line)
                    
def creategraph(folder, graphfile,nthreads, ifhighqual = 1, prior = ""):

    priors = set(prior.split(","))
    
    #allfiles = [x[1] for x in sorted([( 10e20 + os.path.getsize(folder+"/"+file) if  file[:-3] in priors else  5e20 + os.path.getsize(folder+"/"+file) if ( "NC_" in file or "CHM13_" in file) else 10e19 + os.path.getsize(folder+"/"+file) if  ( ("_chr" in file or "HG38_" in file) and 'HG38alt' not in file ) else  10e18 + os.path.getsize(folder+"/"+file) if  ("HG38alt" in file ) else 10e17 + os.path.getsize(folder+"/"+file) if ("HG002_" in file)  else os.path.getsize(folder+"/"+file) , folder+"/"+file ) for file in os.listdir(folder) if file[-3:]==".fa"], reverse = 1)]
   
    allfiles = fileordered([folder + "/"  + file for file in os.listdir(folder) if file.endswith(".fa")], set([x for x in prior.split(",") if len(x)]))

    if len(allfiles) == 0:
            return
    
    editname(allfiles[0],graphfile) 
    
    graphfile_filename = graphfile.split("/")[-1]
    
    dbpath = folder+graphfile_filename+"_db"
    
    if ifhighqual:
        dbcmd = "bash {}/runmakeblastdb -in {} -out {}".format(script_folder,graphfile, dbpath)
    else:
        dbcmd = "ln -f {} {}".format(graphfile, dbpath)
        
    os.system(dbcmd)
    
    for addfile in allfiles[1:]:
        
        insertfiles = addnewseq(folder, graphfile, addfile,nthreads, 0.90, ifhighqual)
        if len(insertfiles):
            os.system(dbcmd)
            
            
def run_align(graphfile, thefile, graphalign, simi, ifhighqual = 1):
    
    lowcasecount = checkmask(thefile)
    
    ifrepeat = 0
    if lowcasecount > 5000:
        ifrepeat = 1
        
    output = thefile + "_align.out"
    
    callalign(thefile,graphfile, output, 1, simi, ifhighqual, ifrepeat)
    queryname, fullpath, fullcigar, pathranges,qranges = graphcigar(graphfile, thefile , output)
    
    lock1.acquire()
    with open(graphalign, mode = 'a') as w:
        w.write("{}\t{}\t{}\t{}\t{}\n".format(queryname, fullpath, fullcigar, pathranges,qranges))
    lock1.release()
    
    
def run_cigar(graphfile, thefile, dbpath):
    
    output = thefile + "_align.out"
    
    queryname, fullpath, fullcigar, pathranges,qranges = graphcigar(graphfile, thefile , output)
    
    return [queryname, fullpath, fullcigar, pathranges,qranges]

def aligngraph(folder, graphfile, graphalign,nthreads,ifblast = 0, ifunfinish = 0):
    
     
    finished = set()
    if os.path.isfile(graphalign):
    
        if ifunfinish:
            iferror = 0
            with open(graphalign, mode = 'r') as f:
                for line in f:
                    line = line.split()
                    if len(line) == 5:
                        name = line[0]
                        finished.add(name+".fa")
                    elif len(line) != 5:
                        iferror = 1
            if iferror:
                with open(graphalign, mode = 'r') as f, open(graphalign+".temp", mode = 'w') as w:
                    for line in f:
                        if len(line.split()) != 5:
                            continue
                        else:
                            w.write(line)
                os.system("mv {} {}".format(graphalign+".temp", graphalign))
        else:
            os.remove(graphalign)

    graphfile_filename = graphfile.split("/")[-1]
    
    dbpath = folder+graphfile_filename+"_db"

    if ifblast:
        dbcmd = "bash {}/runmakeblastdb -in {} -out {} ".format(script_folder,graphfile, dbpath)
        os.system(dbcmd)
    
    dbcmd = "ln -f {} {}".format(graphfile, folder+graphfile_filename)
    os.system(dbcmd)
    
    allfiles = [x[1] for x in sorted([( 10e20 + os.path.getsize(folder+"/"+file) if ( "NC_" in file or "CHM13" in file ) else 10e19 + os.path.getsize(folder+"/"+file) if  ( ("_chr" in file or "HG38" in file) and 'HG38alt' not in file ) else  10e18 + os.path.getsize(folder+"/"+file) if  ("HG38alt" in file )  else os.path.getsize(folder+"/"+file) , folder+"/"+file ) for file in os.listdir(folder) if file.endswith(".fa") and file not in finished], reverse = 1)]
    
    
    with mul.Pool(processes=nthreads) as pool:
        
        pool.starmap(run_align, [(folder+graphfile_filename, thefile, graphalign, 0.9,ifblast) for thefile in allfiles])
        
        
        
        
def main(args):
    
    if len(args.dir)==0:
        folder = args.input +  datetime.datetime.now().strftime("%y%m%d_%H%M%S/")
        #folder = "./tempx/"
        
        folder = folder.replace("%","_")
        folder2 = folder + "/graphtemp/"
        
        try:
            os.system("rm -rf {}".format(folder))
            os.mkdir(folder)
        except:
            pass
            
        try:
            os.system("rm -rf {}".format(folder2))
            os.mkdir(folder2)
        except:
            pass
            
    else:
        folder = args.dir
        folder2 = folder + "/graphtemp/"
        try:
            os.system("rm -rf {}".format(folder2))
            os.mkdir(folder2)
        except:
            pass
            
    if 1 == 1:   
        inputfile_name = args.input.split("/")[-1]
        
        graphfile_raw = folder + inputfile_name+"_graphwithrep.FA"
        graphfile = args.input +"_graph.FA"
        graphalign = args.input + "_allgraphalign.out"
        graphlinear = args.input + "_lineargraph.gaf"
        
        if args.split:
            filesplit(args.input, folder, 0)
            
        if args.create:
            creategraph(folder, graphfile_raw,args.thread, args.ifhighqual, args.prior)
            
        if args.fine:
            runcleanrepeats(folder2, graphfile_raw, graphfile,args.thread, args.ifhighqual)
            
        if args.align:
            aligngraph(folder, graphfile, graphalign,args.thread, args.ifhighqual, args.unfinish)
            
        if args.linear:
            graphtolinear(graphalign, args.input, graphfile, graphlinear, args.thread, folder+"stretcherouts/", args.globalalign)
            
        if args.dir == "":
            pass
            #os.system("rm -rf {} || true".format(folder))
            
            
def run():
    """
            Parse arguments and run
    """
    parser = argparse.ArgumentParser(description="program distract overlap genes and alignments on contigs")
    
    parser.add_argument("-i", "--input", help="path to output file", dest="input", type=str,required=True)
    parser.add_argument("-c", "--create", help="if generate fine graph", dest="create", type=int,default = 1)
    parser.add_argument("-f", "--fine", help="if generate fine graph", dest="fine", type=int,default = 1)
    parser.add_argument("-a", "--align", help="if align to graph", dest="align", type=int,default = 1)
    parser.add_argument("-l", "--linear", help="if align to graph", dest="linear", type=int,default = 1)
    parser.add_argument("-t", "--thread", help="if align to graph", dest="thread", type=int,default = 1)
    parser.add_argument("-s", "--split", help="if align to graph", dest="split", type=int,default = 1)
    parser.add_argument("-d", "--dir", help="if align to graph", dest="dir", type=str,default = "")
    parser.add_argument("-g", "--globalalign", help="if align to graph", dest="globalalign", type=int,default = 0)
    parser.add_argument("-q", "--highqual", help="if align to graph", dest="ifhighqual", type=int,default = 1)
    parser.add_argument("-u", "--unfinish", help="if restart", dest="unfinish", type=int,default = 1)
    parser.add_argument("-p", "--prior", help="prior prefix", dest="prior", type=str,default = "")
    parser.set_defaults(func=main)
    args = parser.parse_args()
    args.func(args)
    
    
if __name__ == "__main__":
    run()
