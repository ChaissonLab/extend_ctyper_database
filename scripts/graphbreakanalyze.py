#!/usr/bin/env python3

import os
import threading
import argparse
import re
import numpy as np
import collections as cl
from scipy import stats


script_folder = os.path.dirname(os.path.abspath(__file__))

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

def readcontigsinfo(inputfile, breakfile):
    
    seqinfo = dict()
    with open(inputfile, mode = 'r') as f:
        for line in f:
            if line.startswith(">"):
                line = line.strip().split()
                name = line[0][1:]
                locus = line[1]
                strd = locus[-1]
                contig, coordi = locus[:-1].split(":")
                start, end = map(int, coordi.split("-"))
                start, end = min(start,end), max(start,end)
                seqinfo[name] = [contig, strd, start, end, line[2], line[3]]
    
    
    locuslocations = dict()
    geneinfo = dict()
    
    with open(breakfile, mode = 'r') as f:
        for index, line in enumerate(f):
            line = line.strip().split()
            contig, strd, start, end = line[:4]
            
            start, end = min(int(start), int(end)), max(int(start), int(end))
            
            locusname = "loci_"+str(index+1)
            
            names = line[4].split(";")
            if len(names) == 1:
                
                
                if abs(max(seqinfo[names[0]][2],seqinfo[names[0]][3]) - end) < 10:
                    geneinfo[contig+"_front"] = [locusname, start, min(seqinfo[names[0]][2],seqinfo[names[0]][3])]
                    print("ledge:",locusname)
                else:
                    geneinfo[contig+"_end"] = [locusname,  max(seqinfo[names[0]][2],seqinfo[names[0]][3]), qend]
                    print("redge:",locusname)
            
            if "NC_0609" in contig:
                locuslocations["loci_"+str(index+1)] = [contig, strd, start, end]
                locusname = "Ref_"+locusname
            
            locuslocations["loci_"+str(index+1)] = [contig, strd, start, end]
            locuslocations[locusname ] = [contig, strd, start, end]
            
            
            
            for name in names:
                if strd == '+':
                    qstart = seqinfo[name][2] - start
                    qend = seqinfo[name][3] - start
                else:
                    qstart = end - seqinfo[name][3]
                    qend = end  - seqinfo[name][2]
                    qstart,qend = min(qstart,qend), max(qstart,qend)
                
                geneinfo[name] = [locusname, qstart, qend]              
    
                
                
    return seqinfo, geneinfo, locuslocations


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
        

        allbreaks =  [break_uniformed.get(x, x) for y in list(oldbreak) for x in y]
        
        allbreaks = sorted( allbreaks + list(addbreaks[locusname]) * 2)
        
        starts = [allbreaks[2*i] for i in range(len(allbreaks)//2)]
        ends = [allbreaks[2*i+1] for i in range(len(allbreaks)//2)]
        
        removelist = set()
        for abreak in list(removebreaks[locusname]):
            
            clean_end_index = 0
            clean_start_index = len(ends)-1
            clean_start_dis = 100000000000
            clean_end_dis = 1000000000000
            for i, start in enumerate(starts):
                if start >= abreak:
                    clean_start_index = i
                    clean_start_dis = start - abreak 
                    break
                    
            for i, end in enumerate(ends):
                if end <= abreak:
                    clean_end_index = i 
                    clean_end_dis = abreak - end 
                else:
                    break
            
            if clean_end_index != clean_start_index -1:
                index = clean_start_index
            else:
                if clean_end_dis > clean_start_dis:
                    index = clean_start_index
                else:
                    index = clean_end_index + 1
                    
            if index not in [0, len(ends)]:
                removelist.add(2*index)
                removelist.add(2*index-1)
            
        allbreaks = [x for i,x in enumerate(allbreaks) if i not in removelist]

        allbreaks_q = gposi_toqposi(allcigars[locusname], allbreaks)
        
        allbreaks_q = sorted(allbreaks_q)
            
        results[locusname] = allbreaks_q
    
    return results
    
    

def filterbreaks(gbreaks, break_contigs, break_count, contigspans):
    
    for name, breaks in contigspans.items():
        
        pass
    

    allcontigs_breaks = [x  for y in gbreaks.values() for x in [ y[0][0]-5, y[-1][1]+5]]
    
    allcontigs_spans = [z  for y in contigspans.values() for x in y for z in x]
    
    allnames =  [z  for name, y in contigspans.items() for x in y for z in [name,name]] 
    
    allbreaks = list(break_count.keys())
    
    break_count_corr = cl.defaultdict(int)
    
    l1, l2 = len(allcontigs_breaks), len(allcontigs_spans) + len(allcontigs_breaks)
    
    allcoordis = allcontigs_breaks + allcontigs_spans+ allbreaks
    
    allcoordis_sort = sorted(list(range(len(allcoordis))), key = lambda x: allcoordis[x])
    
    locus_span = 0
    segment_span = 0
    
    locus_lastposi = cl.defaultdict(int)
    locus_addbreaks = cl.defaultdict(set)
    locus_removebreaks = cl.defaultdict(set)
    allbreaks_found = []
    current_names = set()
    lastbreak = 0
    for sortindex, index in enumerate(allcoordis_sort):
        
        if index < l1:
            if index % 2 == 0:
                locus_span += 1
            else:
                locus_span -= 1
        elif index < l2:
            if index % 2 == 0:
                segment_span += 1
                current_names.add(allnames[index-l1])
            else:
                segment_span -= 1
                current_names.remove(allnames[index-l1])
                
        else:
            coordi = allcoordis[index]
            
            common = break_contigs[coordi] & current_names
            totalbreak = len(common)
            totalcontigs = len(current_names)
            
            ifref = any(x for x in common if x.startswith('Ref_'))
            
            if ifref == 0:
                ifref = -1 if any(x for x in current_names if x.startswith('Ref_')) else 0
            
            
            allbreaks_found.append([coordi,totalbreak, totalcontigs, ifref, current_names ^ common , common]) 
            
    allbreaks_foundsort = sorted(allbreaks_found, key = lambda x: (x[3],x[1]/max(1,x[2])), reverse=1)
    
    existingbreak = [0,max(allcoordis)]
    
    for abreak in allbreaks_foundsort:
        
        coordi, totalbreak, totalcontigs, ifref, uncommon, common = abreak
        
        mindis = min([abs(x-coordi) for x in existingbreak])
        
        ifcut = 0 
        if mindis < 10000:
            ifcut = 0 
        
        elif ifref:
            ifcut = 1 if ifref > 0 else 0
            
        elif mindis > 50000:
            ifcut = 1
            
        elif totalbreak > totalcontigs//2:
            ifcut = 1 
        else:
            ifcut = 0
        
        if ifcut:
            existingbreak.append(coordi)
            for name in uncommon:
                locus_addbreaks[name].add(coordi)
        else:
            for name in common:
                locus_removebreaks[name].add(coordi)
            
            
    return locus_addbreaks, locus_removebreaks


def overlapsegments(oldseg, newseg):
    
    
    
    return overlap
    
    
    

def breaks_ongraph(geneinfo, graphinfo, contigspan):
    
    gbreaks = cl.defaultdict(list)
    
    segmentinfo =  cl.defaultdict(list)
    for name, info in geneinfo.items():
        
        locus,start,end = info
    
        if locus.startswith("Ref_"):
            graphinfo[locus] = graphinfo[locus[4:]]
        
        gstart,gend = qposi_togposi(graphinfo[locus], [start,end-10])
            
        gbreaks[locus].append([gstart,gend+10])
        segmentinfo[locus].append([gstart,gend+10, name])
    
    
    
    gbreaks = {name:sorted(values) for name, values in gbreaks.items()}
    
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
    
    
    return finalbreaks

def annotate_regions(breaksoncontigs, locuslocations, seqinfo):
    
    geneoncontigs = cl.defaultdict(list)
    for genename, info in seqinfo.items():
        
        geneoncontigs[info[0]].append(info[2:]+[genename])
    
    
    regions = []
    for contig, breaks in breaksoncontigs.items():
        
        segments = [[breaks[2*i],breaks[2*i+1]] for i in range(len(breaks)//2)]
        
        oldsegments = geneoncontigs[contig]
        
        for segment in segments:
            
            overlap = sorted([(max(x[1],segment[1]) - min(x[0],segment[0]) - (x[1] - x[0]) - (segment[1] - segment[0]), x) for i,x in enumerate(oldsegments)], reverse = 1)[0]
        
            regions.append([contig, segment[0], segment[1]] + overlap[1][2:])
    
    new_regions = []
    haplo_counter = cl.defaultdict(int)
    for region in regions:
        
        contig, start, end, ifexon, mapinfo,oldname = region
            
        haplo_counter[contig] += 1
        
        newname = "_".join(oldname.split("_")[:-1]) +"_"+ str(haplo_counter[contig])
    
        header = "{}\t{}:{}-{}\t{}\t{}".format(newname, contig, start, end, ifexon, mapinfo  )
    
        new_regions.append([newname,  header])
    
    return new_regions

def locatebreaksoncontigs(breaksonlocus, locuslocations):
    
    locussoncontig = cl.defaultdict(list)
    breaksoncontig = cl.defaultdict(list)
    for locusname,qbreaks in breaksonlocus.items():
            
        
        contig, strd, lstart, lend = locuslocations[locusname]
        
        if strd == "+":
            lbreaks = [lstart + qbreak for qbreak in qbreaks]
        else:
            lbreaks = [lend - qbreak for qbreak in qbreaks]
            
        breaksoncontig[contig].extend([x for x in lbreaks])
        locussoncontig[contig].extend([lstart, lend])
    
    
    results = cl.defaultdict(list)
    breaksoncontig = {contig: sorted(values) for contig, values in breaksoncontig.items()}
    for contig, cbreaks in breaksoncontig.items():
        
        lloci = locussoncontig[contig]
        
        
        themin, themax = min( lloci), max(lloci)
        lastbreak = -10000
        allbreaks = []
        currbreaks = []
        for coordi in cbreaks:
            
            if (coordi - lastbreak) < 1000:
                currbreaks.append(coordi)
            else:
                if len(currbreaks):
                    
                    allbreaks.append(cl.Counter(currbreaks).most_common(1)[0][0])
                currbreaks = [coordi]
            lastbreak = coordi
        
        if len(allbreaks) == 1:
            if  allbreaks[0] - themin > 1000:
                allbreaks.append(themin)
            
            if themax - allbreaks[0] > 1000:
                allbreaks.append(themax)
        
        results[contig] = allbreaks
        
    return results

def uniformbreaks(inputfile, outputfile):
    
    breakfile = outputfile+"_breaks.txt"
    graphref = outputfile+"_breaks.txt.fasta"
    
    linearfile = outputfile +"_breaks.txt.fasta_lineargraph.gaf"
    
    seqinfo, geneinfo, locuslocations = readcontigsinfo(inputfile, breakfile)
    
    pathsize, graphinfo, contigspan= readlineargraph(linearfile)
    
    breaksonlocus = breaks_ongraph(geneinfo, graphinfo, contigspan)
    
    breaksoncontigs = locatebreaksoncontigs(breaksonlocus, locuslocations)
    
    regions = annotate_regions(breaksoncontigs,  locuslocations, seqinfo)
    
    return regions

    
    
            
def main(args):
    
    regions = uniformbreaks(args.input, args.output)
    
    makefasta( regions, args.input, graphref, args.output)
    
    
    
    
    
    
def run():
    """
        Parse arguments and run
    """
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-i", "--input", help="path to input data file",dest="input", type=str, required=True)
    parser.add_argument("-n", "--norm", help="path to input data file",dest="norm", type=str, required="")
    parser.add_argument("-k", "--kmer", help="path to input data file",dest="kmer", type=str, default = "")
    parser.add_argument("-o", "--output", help="path to output file", dest="output",type=str, required=True)
    parser.add_argument("-t", "--threads", help="path to output file", dest="threads",type=int, default = 1)
    parser.add_argument("-q", "--query", help="path to output file", dest="query",type=str, default = "")
    
    parser.set_defaults(func=main)
    args = parser.parse_args()
    args.func(args)
    
    
if __name__ == "__main__":
    run()
