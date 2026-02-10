#!/usr/bin/env python3

import os
import threading
import argparse
import concurrent.futures
import datetime
import collections as cl
import re
import numpy as np


script_folder = os.path.dirname(os.path.abspath(__file__))
cutoffdistance = 30000

def process_locus2(region, haplopath, folder, outputfile):
    
    newname,  contig, start, end, ifexon, mapinfo = region
    
    if "#" not in contig: 
        haplo = "CHM13_h1" if "NC_0609" in contig else "HG38_h1"
    else:
        haplo = "_h".join(contig.split("#")[:2])
        
    path = haplopath[haplo]
    
    outputfile0 = os.path.join(folder, f"{newname}.fa")
    
    header = ">{}\t{}:{}-{}\t{}\t{}".format(newname, contig, start, end, ifexon, mapinfo)
    
    cmd1 = f'echo "{header}" > {outputfile0} && samtools faidx {path} {contig}:{start}-{end} | tail -n +2 >> {outputfile0}'
 
    os.system(cmd1)
    
    return outputfile0


def makenewfasta2(regions, queryfile, outputfile, folder, threads):
    
    haplopath = dict()
    with open(queryfile, mode = 'r') as f:
        for line in f:
            line = line.strip().split()
            haplopath[line[0]] = line[1]
            
    os.system("echo > "+outputfile)
    
    alloutputs = []

    if 1==1:
        with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
            futures = []
            for region in regions:
                # Submit jobs to executor
                #process_locus2( region, haplopath, folder, outputfile)
                futures.append(executor.submit(process_locus2, region, haplopath, folder, outputfile))
                alloutputs.append(os.path.join(folder, f"{region[0]}.fa"))
            # Wait for all threads to complete
            #concurrent.futures.wait(futures)
            
    for result in alloutputs:
        
        outputfile0 = result
        cmd2 = f"cat {outputfile0} >> {outputfile}"
        os.system(cmd2)
        

def makegraph(fastafile, folder, threads):
    print("python {}/graphmake.py -i {}  -d {} -t {} -l 0 -c 0 -f 0 ".format(script_folder, fastafile, folder, threads))
    os.system("python {}/graphmake.py -i {}  -d {} -t {} -l 0 -c 0 -f 0 ".format(script_folder, fastafile, folder, threads))
    
def process_locus(index, line, haplopath, folder, outputfile):
    line = line.strip().split()
    contig = line[0]
    
    if "#" not in contig: 
        haplo = "CHM13_h1" if "NC_0609" in contig else "HG38_h1"
    else:
        haplo = "_h".join(contig.split("#")[:2])
        
    path = haplopath[haplo]
    strd = "-i" if line[1] == "-" else ""
    
    outputfile0 = os.path.join(folder, f"loci_{index}.fa")
    
     
    header = ">loci_{}\t{}:{}-{}{}\t{}".format(index, contig, line[2], line[3], line[1], line[4])
    
    line[2] = str(int(line[2]) + 1) 
    cmd1 = f'echo "{header}" > {outputfile0} && samtools faidx {path} {contig}:{line[2]}-{line[3]} {strd} | tail -n +2 >> {outputfile0}'
    os.system(cmd1)
    
    return outputfile0

    
def makenewfasta(locifile, queryfile, outputfile, folder, threads):
    
    haplopath = dict()
    with open(queryfile, mode = 'r') as f:
        for line in f:
            line = line.strip().split()
            haplopath[line[0]] = line[1]
            
    os.system("echo > "+outputfile)
    
    alloutputs = []
    with open(locifile, mode='r') as f:
        with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
            futures = []
            for index, line in enumerate(f, start=1):
                # Submit jobs to executor
                #process_locus(index, line, haplopath, folder, outputfile)
                futures.append(executor.submit(process_locus, index, line, haplopath, folder, outputfile))
                alloutputs.append(os.path.join(folder, f"loci_{index}.fa"))
            # Wait for all threads to complete
            concurrent.futures.wait(futures)
            
    for result in alloutputs:
        
        outputfile0 = result
        cmd2 = f"cat {outputfile0} >> {outputfile}"
        os.system(cmd2)
        
        
def findbreaks(inputfile, outputfile):
    
    ifexons = set()
    names = []
    contigs = cl.defaultdict(list)
    name_tocontig = cl.defaultdict(list)
    index = -1
    with open(inputfile, mode = 'r') as f:
        
        for line in f:
            
            if line.startswith(">"):
                
                index += 1
                name, locus = line[1:].split()[:2]
                names.append(name)
                exon = line.split()[-2]
                if exon == "Exon":
                    ifexons.add(name)
                    
                strd = "+"
                contig, region = locus.split(":")
                if region[-1] in ["+","-"]:
                    region = region[:-1]
                    strd = locus[-1]
                    
                region = region.split("-")
                region = [int(region[0]), int(region[1]), name] if strd == "+" else [int(region[1]), int(region[0]), name]
                
                contigs[contig].append(region)
                
    allsize =0 
    new_contigs = cl.defaultdict(list)
    for contig, regions in contigs.items():
        
        new_regions = []
        
        strd = sorted(regions, key = lambda x:  -abs(x[0]-x[1]))[0]
        
        strd = '+' if strd[1] > strd[0] else '-'
        
        regions_sort = sorted(regions) if contig[-1] == '+' else sorted(regions, reverse = 1)
        
        lastend = regions_sort[0]
        lastregion = regions_sort[0]
        new_regions = [[regions_sort[0]]]
        
        allsize += abs(regions_sort[0][1]-regions_sort[0][1])
        for region in regions_sort[1:]:
            if max(region[0],region[1], lastend[0],lastend[1]) - min(region[0],region[1], lastend[0],lastend[1]) - abs(region[1]-region[0] ) - abs(lastend[1] - lastend[0]) < 30000:
                
                new_regions[-1].append(region)
            else:
                new_regions.append([region])
            
            allsize += abs(region[0]-region[1])
            lastend = region
                
        new_contigs.update({(contig+strd,i):regions for i,regions in enumerate(new_regions)})
    

    with open(outputfile, mode = 'w') as w:
        for (scarf, index), regions in new_contigs.items():
            locus = "{}".format(scarf[:-1])
            start = min(regions[0][0],regions[0][1],regions[-1][0],regions[-1][1])
            end = max(regions[0][0],regions[0][1],regions[-1][0],regions[-1][1])
            names = ";".join([region[-1] for region in regions])

            w.write("{}\t{}\t{}\t{}\t{}\n".format(locus,scarf[-1],max(0,start-30000),(end+30000),names))
            
            
    return new_contigs


def findstartend(cigar_str):
    
    if not cigar_str:
        return 0,0
    
    cigars = re.findall(r'\d+[M=XIHD]', cigar_str)
    
    lastposi=0
    lsize = 0
    rsize = 0
    for i,cigar in enumerate(cigars):
        
        curr_size = int(cigar[:-1])
        char = cigar[-1]
        
        if (char == '=' or char ==  'M' or char == 'I' ) and curr_size>50:
            break
        if char in ['D','H','X','M','=']:
            lsize += curr_size
            
    for i,cigar in enumerate(cigars[::-1]):
        
        curr_size = int(cigar[:-1])
        char = cigar[-1]
        
        if (char == '=' or char ==  'M' or char == 'I' ) and curr_size>50:
            break
        if char in ['D','H','X','M','=']:
            rsize += curr_size
            
    return lsize,rsize


def qposi_togposi(cigar_str, find_qposes):
    
    l = len(find_qposes)
    if l == 0:
        
        return []
    
    allcigars = re.findall(r'[><][^><]+', cigar_str)
    
    
    curr_rposi = 0
    curr_qposi = 0
    new_rposi = 0
    new_qposi = 0
    
    curr_strd = 1
    curr_size = 0
    
    find_qposes = find_qposes + [0]
    find_qpos_index = 0
    find_qpos = find_qposes[0]
    find_rposes = []
    
    for cigars in allcigars:
        
        curr_strd = 1 if cigars [0]=='>' else -1
        
        curr_posi, cigar  = cigars.split(":")
        
        curr_rposi,curr_qposi = curr_posi[1:].split("_")
        
        curr_rposi,curr_qposi = int(curr_rposi), int(curr_qposi)
        
        cigars = re.findall(r'\d+[M=XIHD]', cigars)
        
        new_rposi, new_qposi = curr_rposi,curr_qposi
        
        for cigar in cigars:
            
            curr_size = int(cigar[:-1])
            char = cigar[-1]
            
            if char == '=' or char ==  'M':
                
                new_rposi = curr_rposi + curr_strd * curr_size
                new_qposi = curr_qposi + curr_size
                
            elif char == 'I':
                
                new_qposi = curr_qposi + curr_size 
                
            elif char == 'D' or char == 'H':
                
                new_rposi = curr_rposi + curr_strd * curr_size
                
            elif char == 'X':
                
                new_rposi = curr_rposi + curr_strd * curr_size
                new_qposi = curr_qposi + curr_size
                
            elif char == 'S':
                
                if curr_qposi == 0:
                    curr_qposi = curr_size
                    
                curr_size = 0
                continue
            
            while find_qpos_index<l and find_qpos <= new_qposi:
                
                
                if char == '=' or char ==  'M' or char =='X' :
                    
                    find_rposes.append(curr_rposi + curr_strd * (find_qpos-curr_qposi))
                else:
                    
                    find_rposes.append(curr_rposi)
                    
                find_qpos_index += 1
                find_qpos = find_qposes[find_qpos_index]
                
            if find_qpos_index == l:
                return find_rposes
            
            curr_size = 0
            curr_rposi = new_rposi
            curr_qposi = new_qposi
            
            
    while find_qpos_index<l:
        
        find_qpos_index += 1
        find_rposes.append(curr_rposi)
        
    return find_rposes


def gposi_toqposi(cigar_str, find_rposes):
    
    l = len(find_rposes)
    if l == 0:
        
        return []
    
    allcigars = re.findall(r'[><][^><]+', cigar_str)
    allcigars = sorted(allcigars, key = lambda x: int(x.split(":")[0].split("_")[0][1:]))
    
    
    find_rposes = find_rposes + [0]
    curr_rposi = 0
    curr_qposi = 0
    new_rposi = 0
    new_qposi = 0
    
    curr_size = 0
    curr_strd = 1
    find_rpos_index = 0
    find_rpos = find_rposes[0]
    find_qposes = []
    
    for cigars in allcigars:
        
        new_strd = 1 if cigars[0]=='>' else -1
        new_posi, cigars  = cigars.split(":")
        
        new_rposi,new_qposi = new_posi[1:].split("_")
        
        new_rposi,new_qposi = int(new_rposi), int(new_qposi)
        
        rsize = sum([int(x[:-1]) for x in re.findall(r'\d+[M=XHD]', cigars)]+[0])
        qsize = sum([int(x[:-1]) for x in re.findall(r'\d+[M=XI]', cigars)]+[0]) 
        
        if new_strd == 1:
            cigars = re.findall(r'\d+[M=XHDI]', cigars)
            
        else:
            cigars = re.findall(r'\d+[M=XHDI]', cigars)[::-1]
            new_rposi -= sum([int(x[:-1]) for x in cigars if x[-1] in "M=XHD"]+[0])
            new_qposi += sum([int(x[:-1]) for x in cigars if x[-1] in "M=XI"]+[0])
            
            
            
        while find_rpos_index<l and find_rpos <= new_rposi + (1-curr_strd)-1//2 :
            
            if abs(find_rpos - curr_rposi) < abs(find_rpos - new_rposi):
                find_qposes.append(curr_qposi)
            else:
                find_qposes.append(new_qposi)
                
            find_rpos_index += 1
            find_rpos = find_rposes[find_rpos_index]
            
        if find_rpos_index == l:
            return find_qposes
        
        curr_strd, curr_rposi, curr_qposi = new_strd, new_rposi, new_qposi
        
        
        for cigar in cigars:
            
            curr_size = int(cigar[:-1])
            char = cigar[-1]
            
            if char == '=' or char ==  'M':
                
                new_rposi = curr_rposi + curr_size 
                new_qposi = curr_qposi + curr_size * curr_strd
            elif char == 'I':
                
                new_qposi = curr_qposi + curr_size * curr_strd
                
            elif char == 'D' or char == 'H':
                
                new_rposi = curr_rposi + curr_size 
                
            elif char == 'X':
                
                new_rposi = curr_rposi + curr_size 
                new_qposi = curr_qposi + curr_size * curr_strd
                
                
            elif char == 'S':
                
                if curr_qposi == 0:
                    curr_qposi = curr_size 
                    
                curr_size = 0
                continue
            
            while find_rpos_index<l and find_rpos <= new_rposi + (1-curr_strd)-1//2 :
                
                
                if char == '=' or char ==  'M' or char =='X' :
                    
                    find_qposes.append(curr_qposi + curr_strd * (find_rpos-curr_rposi))
                else:
                    
                    find_qposes.append(curr_qposi)
                    
                find_rpos_index += 1
                find_rpos = find_rposes[find_rpos_index]
                
            if find_rpos_index == l:
                return find_qposes
            
            curr_size = 0
            curr_rposi = new_rposi
            curr_qposi = new_qposi
            
            
    while find_rpos_index<l:
        
        find_rpos_index += 1
        find_qposes.append(curr_qposi)
        
    return find_qposes



def getspan_ongraph(path, cigar):
    
    span = []
    cigar = cigar.translate(str.maketrans('', '', 'ATCGatcgNn-'))
    cigars = re.findall(r'[><][^><]+', cigar)
    pathes = re.findall(r'[><][^<>]+',path)
    
    
    newcigars = cl.defaultdict(str)
    lcut, rcut  =0 ,0
    rposi = 0
    for cigarstr in cigars:
        
        thename,thecigar = cigarstr.split(":")
        
        asize = sum([int(x[:-1]) for x in re.findall(r'\d+[M=]', thecigar)]+[0])
        
        rsize = sum([int(x[:-1]) for x in re.findall(r'\d+[M=XHD]', thecigar)]+[0])
        
        qsize = sum([int(x[:-1]) for x in re.findall(r'\d+[M=XI]', thecigar)]+[0])
        
        newcigars[thename[1:]] = (thecigar, rposi, asize, rsize, qsize)
        
        if asize >= 1000:
            
            lcut, rcut = findstartend(thecigar)
            if lcut+rcut > rsize:
                continue

            if len(span) and rposi + lcut - span[-1][1] < 2000:
                span[-1][1] = rposi+rsize - rcut
            else:
                span.append([rposi+ lcut,rposi+rsize- rcut])
                
        rposi += rsize
        
        
    qposi = 0
    cigar_order = []
    reserved = ""
    reserved_size = 0
    lastaligned = -1
    for apath in pathes:
        
        thecigar, rposi, asize, rsize, qsize = newcigars[apath[1:]]
        
        if asize < 1000:
            
            thecigar = "{}H".format(rsize)
            
            if qsize >0:                
                if lastaligned >= 0:
                    cigar_order[lastaligned] += "{}I".format(qsize)
                    qposi += qsize
                else:
                    reserved = "{}I".format(qsize)
                    reserved_size = qsize
        else:
            if apath.startswith("<"):
                
                thecigar = "".join(re.findall(r'\d+[M=XHDI]', thecigar)[::-1])
                rposi += rsize
                
            if len(reserved):
                thecigar = reserved + thecigar
                reserved = ""
                qsize += reserved_size
                reserved_size = 0
                
            cigar_order.append(apath[0]+str(rposi)+"_"+str(qposi )+":"+thecigar)
            lastaligned = len(cigar_order) - 1
            
            qposi += qsize 
            
            
    newcigars= "".join(cigar_order)
    
    return "".join(newcigars), span



def combinebreaks(gbreaks):
    
    breaks = [x for y in gbreaks.values() for x in y]
    names = [name for name,y in gbreaks.items() for x in y]
    
    all_coordinates = [x for y in breaks for x in y[:2]]
    
    sort_index = sorted(range(len(all_coordinates)), key = lambda x: all_coordinates[x])
    
    break_groups = []
    break_group = []
    
    break_group_index = set([])
    break_group_indexs = []
    
    last_coordinate = 0
    for index in sort_index:
        
        coordinate = all_coordinates[index]
        
        if coordinate - last_coordinate < 300:
            
            break_group.append(coordinate)
            break_group_index.add(names[index//2])
            
        else:
            
            break_groups.append(break_group)
            break_group_indexs.append(break_group_index)
            
            break_group = [coordinate]
            break_group_index = set([names[index//2]])
            
        last_coordinate = coordinate
        
    break_groups.append(break_group)
    break_group_indexs.append(break_group_index)
    
    extendgroup = cl.defaultdict(set)
    
    
    return break_groups, break_group_indexs, extendgroup


def readlineargraph(linearfile):
    
    contigspan = dict()
    pathsize = dict()
    graphinfo = dict()
    with open(linearfile, mode = 'r') as f:
        
        for line in f:
            
            if line.startswith("S") :
                line = line.split()
                name, size = "_".join(line[1].split("_")[1:]), int(line[4].split(":")[-1])
                pathsize[name] = size
                
            else:
                line = line.strip().split()
                
                newcigar, span = getspan_ongraph(line[3],line[5])
 
                graphinfo[line[1]] = newcigar
                
                contigspan[line[1]] = span
                
    return pathsize, graphinfo, contigspan



def applybreaks(oldbreaks, addbreaks, removebreaks, allcigars, break_uniformed):
    
    results = dict()
    for locusname, oldbreak in oldbreaks.items():

        removerange = [ x for y in removebreaks[locusname] for x in [y-1000, y+1000]]
        
        breaks =  sorted([break_uniformed.get(x, x) for y in list(oldbreak) for x in y])
        
        lremove = len(removerange)
        allbreaks = removerange + breaks
        
        allbreaks_sort = [(i,value) for i, value in sorted(enumerate(allbreaks), key=lambda x: x[1])]
       

        lastcoordi = -100000
        removerange_cover = 0
        removelist = set()
        for sortindex,(index,coordi) in enumerate(allbreaks_sort):
            
            if index < lremove:
                if index % 2:
                    removerange_cover -= 1
                else:
                    removerange_cover += 1
            else:
                distance = coordi - lastcoordi 
                if sortindex and distance > 100000:
                    continue
                if sortindex and distance < 10000:
                    removelist.add(index)
                else:
                    if removerange_cover > 0:
                        removelist.add(index)
            
            lastcoordi = coordi
                
        newbreak = [x for i,x in enumerate(breaks) if i not in removelist]
        
        if len(newbreak) <= 1:
            newbreak = [x for y in oldbreak for x in y]

        
        breaks_q = gposi_toqposi(allcigars[locusname], newbreak)
       
        breaks_q = sorted(breaks_q)
        
        results[locusname] = breaks_q
        
    return results



    

def overlapsegments(oldseg, newseg):
    
    
    return overlap


    
    
def breaks_ongraph(genetolocus, graphinfo, contigspan):
    
    gbreaks = cl.defaultdict(list)
    
    segmentinfo =  cl.defaultdict(list)
    for name, infos in genetolocus.items():
        
        for info in infos:

            locus,start,end = info
      
            if locus.startswith("Ref_"):
                graphinfo[locus] = graphinfo[locus[4:]]
           
            gstart,gend = qposi_togposi(graphinfo[locus], [start,end-10])
       
            gbreaks[locus].append([gstart,gend+10])
            segmentinfo[locus].append([gstart,gend+10, name])

    gbreaks = {name:sorted([sorted(x) for x in values]) for name, values in gbreaks.items()}
   
    break_groups, break_group_names, extendgroup = combinebreaks(gbreaks)
    
    break_count = cl.defaultdict(int)
    break_uniformed = dict()
    break_contigs = cl.defaultdict(set)
    for group,names in zip(break_groups,break_group_names):
        
        if len(group) == 0:
            continue
        themedian = int(np.median(group))
        break_contigs[themedian] = names
        break_count[themedian] = len(extendgroup[themedian])
        #print(themedian, len(extendgroup[themedian]), len(group))
        for coordi in group:
            
            break_uniformed[coordi] = themedian
            
            
    gbreaks = {locus: sorted( [[break_uniformed[x[0]], break_uniformed[x[1]]] for x in y ]) for locus,y in gbreaks.items()}
    
    addbreaks, removebreaks = filterbreaks(gbreaks, break_contigs, break_count, contigspan)
    
    finalbreaks = applybreaks(gbreaks, addbreaks, removebreaks, graphinfo, break_uniformed)
    
    
    return gbreaks, finalbreaks

def annotate_regions(breaksoncontigs, locuslocations, genetoscaff):
    
    geneoncontigs = cl.defaultdict(list)
    for genename, info in genetoscaff.items():
        geneoncontigs[info[0]].append(info[2:]+[genename])
        
    regions = []
    for contig, breaks in breaksoncontigs.items():
        
        segments = [[breaks[i],breaks[i+1]] for i in range(len(breaks) - 1)]
        
        oldsegments = geneoncontigs[contig]
        for segment in segments:

            overlap = sorted([(max(x[1],segment[1]) - min(x[0],segment[0]) - (x[1] - x[0]) - (segment[1] - segment[0]), x) for i,x in enumerate(oldsegments)], reverse = 1)[0]
            regions.append(["_".join(contig.split("_")[:-1]), segment[0], segment[1]] + overlap[1][2:])
    
            
    new_regions = []
    haplo_counter = cl.defaultdict(int)
    for region in regions:
        
        contig, start, end, ifexon, mapinfo,oldname = region
        
        haplo_counter[contig] += 1
        
        newname = "_".join(oldname.split("_")[:-1]) +"_"+ str(haplo_counter[contig])
        
        #header = "{}\t{}:{}-{}\t{}\t{}".format(newname, contig, start, end, ifexon, mapinfo  )
        
        new_regions.append([newname,  contig, start, end, ifexon, mapinfo])
        
    return new_regions

def coordinate_uniform(cbreaks_sort):

    uniform_coordis = []
    
    currbreaks =  []
    lastbreak = -1000
    for coordi in cbreaks_sort:
        if (coordi - lastbreak) < 1000:
            currbreaks.append(coordi)
        else:
            if len(currbreaks):
                uniform = cl.Counter(currbreaks).most_common(1)[0][0]
                uniform_coordis.extend([uniform] * 1)
                
            currbreaks = [coordi]
        lastbreak = coordi
    
    if len(currbreaks):
        uniform = cl.Counter(currbreaks).most_common(1)[0][0]  
        uniform_coordis.extend([uniform] * 1)

    return uniform_coordis

def locatebreaksoncontigs(breaksonlocus, locuslocations, oldbreaks):
    
    locusnamecontig = cl.defaultdict(list)
    locussoncontig = cl.defaultdict(list)
    breaksoncontig = cl.defaultdict(list)
    for locusname,qbreaks in breaksonlocus.items():
       
        contig, strd, lstart, lend = locuslocations[locusname]
        
        if strd == "+":
            lbreaks = [lstart + qbreak for qbreak in qbreaks]
        else:
            lbreaks = [lend - qbreak for qbreak in qbreaks]
        
        breaksoncontig[contig].extend(lbreaks)
        locussoncontig[contig].extend([lstart, lend])
        locusnamecontig[contig].extend(locusname)
        
    results = cl.defaultdict(list)
    for scafford, cbreaks in breaksoncontig.items():
                
        cbreaks_sort_index = sorted(range(len(cbreaks)), key = lambda x: cbreaks[x])
        
        cbreaks_sort = [cbreaks[x] for x in cbreaks_sort_index]
        
        cbreaks_sort_uniform = coordinate_uniform(cbreaks_sort)
            
        results[scafford] = cbreaks_sort_uniform 
    
    return results

def uniformbreaks(inputfile, outputfile):
    
    breakfile = outputfile+"_loci.txt"
    graphref = outputfile+"_loci.txt.fasta"
    
    linearfile = outputfile +"_loci.txt.fasta_lineargraph.gaf"
    
    genetoscaff, genetolocus, locuslocations = readcontigsinfo(inputfile, graphref)
   
    pathsize, graphinfo, contigspan= readlineargraph(linearfile)

    oldbreaks, newbreaksonlocus = breaks_ongraph(genetolocus, graphinfo, contigspan)
    
    breaksoncontigs = locatebreaksoncontigs(newbreaksonlocus, locuslocations, oldbreaks)
  
    regions = annotate_regions(breaksoncontigs,  locuslocations, genetoscaff)
 
    return regions
    
def main(args):
    
    #os.system("python {}/exambreaks.py -i {} -n {} -o {} ".format(script_folder,args.input, args.norm, args.input))
    
    if len(args.folder)==0: 
        folder = args.input +  datetime.datetime.now().strftime("%y%m%d_%H%M%S/").replace("%","_")
        os.system("rm -rf {} || true".format(folder))
        os.system("mkdir {} || true ".format(folder))
    else:
        folder = args.folder
        os.system("mkdir {} || true".format(folder))
        
    graphref = args.output+"_loci.txt.fasta"
    
    #pairs = findbreaks(args.input, args.output+"_loci.txt")

    ###
    #makenewfasta(args.output+"_loci.txt", args.query, args.output+"_loci.txt.fasta", folder ,args.threads)
    
    if len(args.kmer):
        pass
        #os.system("python {}/kmerstrd.py -i {} -k {} -o {}".format(script_folder,graphref, args.kmer, graphref+"_"))
        #os.system("mv {} {} || true".format(graphref+"_", graphref))
        
    makegraph(args.output+"_loci.txt.fasta", folder, args.threads)
    
    #regions = uniformbreaks(args.input, args.output)
   
    
    #makenewfasta2(regions, args.query, args.output, folder, args.threads)

    if len(args.folder)==0:
        os.system("rm  -rf {} || true".format(folder))

    """
    if len(args.kmer):
        
        os.system("python {}/kmerstrd.py -i {} -k {} -o {}".format(script_folder,args.output, args.kmer, args.output+"_"))
        os.system("mv {} {} || true".format(args.output+"_", args.output))
    
    if len(args.folder)==0:
        os.system("rm  -rf {} || true".format(folder))
    """
    
    
        
def run():
    """
        Parse arguments and run
    """
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-i", "--input", help="path to input data file",dest="input", type=str, required=True)
    parser.add_argument("-n", "--norm", help="path to input data file",dest="norm", type=str, required="")
    parser.add_argument("-k", "--kmer", help="path to input data file",dest="kmer", type=str, default="")
    parser.add_argument("-o", "--output", help="path to output file", dest="output",type=str, required=True)
    parser.add_argument("-t", "--threads", help="path to output file", dest="threads",type=int, default = 1)
    parser.add_argument("-q", "--query", help="path to output file", dest="query",type=str, default = True)
    parser.add_argument("-d", "--folder", help="path to output file", dest="folder",type=str, default = "")   
    
    parser.set_defaults(func=main)
    args = parser.parse_args()
    args.func(args)
    
    
if __name__ == "__main__":
    run()
