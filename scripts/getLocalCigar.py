#!/usr/bin/env python3

import re
import argparse

regex = r'(\d+)([HSMIDX=])([ATCGNatcgn]*)'

def cigar_cutqrange(fullcigar, qranges):    
    
    def finish_currpath(path, cigars, rsize, strand = 1):
        
        if len(cigars) == 0:
            return ""
        
        suffix = f"{rsize}H" if rsize else ""
        
        if strand == -1:
            cigars = cigars[::-1]
            path = path.replace('<','`').replace('>','<').replace('`','>')
            
        return path+":"+"".join(cigars)+suffix
    
    
    #distract graphic cigar using query coordinates using base-0, right open
    
    allcigars = re.findall(r'[><][^><]+', fullcigar)
    
    total_q = len(qranges)
    
    if total_q == 0:
        return []
    
    results = [[] for x in range(total_q)]
    results_strd = [x[2] if len(x) > 2 else 1  for x in qranges]
    
    curr_qindex = 0
    curr_qregion = qranges[0]
    curr_results = results[0]
    qstart,qend,qstrand = curr_qregion[0],curr_qregion[1],curr_qregion[2] if len(curr_qregion) > 2 else 1
    
    qposi = 0
    
    
    for cigars in allcigars:
        
        strd = cigars[0]
        
        path, cigars = cigars.split(":")
        
        
        cigars = [(int(L), O, S) for L, O, S in re.findall(regex, cigars)]
        rsize = sum([x[0] for x in cigars if x[1] in ['M','X','H','D','=']])
        rposi = 0
        
        newcigars = []
        for cigar in cigars:
            
            lastrposi = rposi
            lastqposi = qposi
            
            thesize,thetype,theseq  = cigar
            
            if thetype == "=" or thetype == "M":
                
                qposi += thesize
                rposi += thesize
                
            elif thetype == "X" or thetype == "S":
                
                qposi += thesize
                rposi += thesize
                
            elif thetype == "D" or thetype == "H":
                
                rposi += thesize
                
            elif thetype == "I":
                
                qposi += thesize
                
                
            while qstart < qposi :
                
                if len(newcigars) == 0 and thetype != 'H':
                    
                    if thetype != 'I' and qstart > lastqposi:
                        newrposi = lastrposi+qstart-lastqposi
                    else:
                        newrposi = lastrposi
                    
                    if newrposi:
                        newcigars.append(f"{newrposi}H")
                    
                segmentsize =  min(qposi, qend)- max(lastqposi, qstart)
                if segmentsize >= 0:
                    if thetype == 'I' and segmentsize < thesize:
                        newseq = f"{theseq[:2]}{theseq[(max(0,qstart-lastqposi)+2):(qend-lastqposi+2)]}"
                        newcigars.append(f"{segmentsize}I{newseq}")
                    elif thetype == "D" or thetype == "H":
                        newcigars.append("{}{}{}".format(thesize,thetype,theseq))
                    else:
                        newcigars.append("{}{}{}".format(segmentsize,thetype,theseq))
                        
                if qposi >= qend:
                    curr_results.append(finish_currpath(path,newcigars,rsize-rposi,qstrand))
                    
                    curr_qindex += 1
                    if curr_qindex >= total_q:
                        
                        for i in range(0,total_q):
                            
                            if results_strd[i] == -1:
                                results[i] = "".join(results[i][::-1])
                            else:
                                results[i] = "".join(results[i])
                                
                        return results
                    
                    newcigars = []
                    curr_qregion = qranges[curr_qindex]
                    curr_results = results[curr_qindex]
                    qstart,qend,qstrand = curr_qregion[0],curr_qregion[1],curr_qregion[2] if len(curr_qregion) > 2 else 1
                else:
                    break
                
        if len(newcigars):
            curr_results.append(finish_currpath(path,newcigars,0,qstrand))
            
    for i in range(0,total_q):
        if results_strd[i] == -1:
            results[i] = "".join(results[i][::-1])
        else:
            results[i] = "".join(results[i])
            
    return  results


def parse_coord_token(token: str):
    
    """
    Parse a coord token like:
        HG02583#1#CM089492.1:85255176-85426487+
    or  HG02583#1#CM089492.1:85255176-85426487
    Returns: (ref_tag, start, end, strand) with strand = +1 / -1
    """
    strand = 1
    if token[-1] in "+-":
        strand = 1 if token[-1] == "+" else -1
        core = token[:-1]
    else:
        core = token
        
    if ":" not in core or "-" not in core:
        raise ValueError(f"Cannot parse coord token: {token}")
        
    left, rng = core.split(":", 1)
    start_s, end_s = rng.split("-", 1)
    return left, int(start_s), int(end_s), strand

def getoverlaps(cols, input_regions):
    
    # merge header like:
    #   >loci_1 HG02583#1#CM089492.1:85255176-85426487+ group50231p1_HG02583_h1_1
    merge_names = cols[0]
    merge_coord_token = cols[1]
    
    # coordinates for the merge region
    contig, m_start, m_end, m_strand = parse_coord_token(merge_coord_token)
    
    overlap_names = [x for x,reg in input_regions.items() if reg["contig"] == contig and max(m_end, reg["end"]) - min(m_start,reg["start"]) + 10 < (reg["end"]-reg["start"]) + (m_end - m_start)]
    
    # full graphic cigar for this merge region
    
    outputs = []
    for in_name in overlap_names:
        if in_name not in input_regions:
            # header says it overlaps, but we didn't see this input name
            continue
        
        reg = input_regions[in_name]
        in_start = reg["start"]
        in_end   = reg["end"]
        in_strand = reg["strand"]
        
        # 6. local query coords: input region relative to merge region
        #    (assumes same coordinate system; adjust if you use 1-based vs 0-based)
        q0 = in_start - m_start
        q1 = in_end   - m_start
        
        # If no overlap at all, skip
        merge_len = m_end - m_start
        if q1 <= 0 or q0 >= merge_len:
            continue
        
        # clamp to merge region if slightly out of bounds
        q0 = max(0, q0)
        q1 = min(merge_len, q1)
        
        # strand for cigar_cutqrange:
        # if input and merge strands differ, use -1
        strand_param = 1 if in_strand == m_strand else -1
    
        outputs.append([q0,q1,strand_param,in_name,  (reg["end"]-reg["start"]) + (m_end - m_start) - max(m_end, reg["end"]) + min(m_start,reg["start"]) ])
    
    return outputs

def getpathinfo(sub_full):

    # 7. Build output columns
    # (b & c) split graphic cigar into paths and cigars
    segments = re.findall(r'[><][^><]+', sub_full)
    if not segments:
        return []
    
    path_list  = []
    cigar_list = []
    r_spans    = []
    q_spans    = []
    
    qposi = 0
    for seg in segments:
        # seg looks like ">path:...cigar..."
        path, cig = seg.split(":", 1)
        path_list.append(path)
        cigar_list.append(cig)
        
        # (d) compute r_span for this path
        ops = [(int(L), O, S) for L, O, S in re.findall(regex, cig)]
    
        rsize = sum(L for L, O, S in ops if O in ["M", "X", "H", "D", "="])
        qsize = sum(L for L, O, S in ops if O in ["M", "X", "I", "="])
        first = ops[0]
        start_r = first[0] if first[1] == "H" else 0
        end_r = start_r + rsize
        r_spans.append(f"{start_r}_{end_r}")
            
        # (e) q_span: (qstart, qend) for this cigar.
        # Here we use absolute input coords from the input header.
        q_spans.append(f"{qposi}_{qposi+qsize}")
        qposi += qsize
        
    paths_str  = "".join(path_list)
    cigar_str  = "".join(cigar_list)
    r_span_str = ";".join(r_spans)
    q_span_str = ";".join(q_spans)
    
    return paths_str,cigar_str,r_span_str,q_span_str

def readinput(inputfile):
    input_names = []
    input_regions = {}
    with open(inputfile) as f:
        for ln in f:
            if not ln.startswith(">"):
                continue
            header = ln[1:].strip()
            cols = header.split()
            
            # name
            name = cols[0]
            
            # coord token: prefer 3rd column, otherwise first token with ":" and "-"
            coord_token = cols[1]
            
            ref_tag, start, end, strand = parse_coord_token(coord_token)
            input_names.append(name)
            input_regions[name] = {
                "ref": ref_tag,
                "start": start,
                "end": end,
                "strand": strand,
                "contig": ref_tag
            }

    return input_names,input_regions

def readaligns(alignfile):
    
    align = {}
    with open(alignfile) as f:
        for ln in f:
            ln = ln.strip()
            if not ln or ln.startswith("#"):
                continue
            parts = ln.split()
            if len(parts) < 3:
                continue
            name, cigar = parts[0], parts[2]
            align[name] = cigar
            
    return align


def main(args):
    # 4. Read align file: first col = merge name, second col = full graphic cigar
    align = readaligns(args.align)
    input_names,input_regions = readinput(args.input)
    
    # 5–7. Read merge FASTA headers, look up their cigars,
    #      compute local query ranges for each overlapping input, call cigar_cutqrange,
    #      and write the new alignment-like output.
    
    overlapsize = dict()
    outputcigars = dict()
    with open(args.merge) as f:
        for ln in f:
            if not ln.startswith(">"):
                continue
            header = ln[1:].strip()
            cols = header.split()
            if len(cols) < 3:
                continue
            
            overlaps = sorted(getoverlaps(cols, input_regions))
            
            overlaps_cigar = cigar_cutqrange(align[cols[0]],overlaps)
            
            for overlap,cigar in zip(overlaps,overlaps_cigar):
            
                if overlap[3] in overlapsize and overlap[4] <= overlapsize[overlap[3]]:
                    continue
                
                overlapsize[overlap[3]] = overlap[4]
                outputcigars[overlap[3]] = cigar
    
    with open(args.output, "w") as out:
        # output header
        #out.write("\t".join(["name", "paths", "cigar", "r_span", "q_span"]) + "\n")
        
        for name in input_names:
            
            cigar_str = outputcigars[name]
            
            paths_str,cigar,r_span_str,q_span_str = getpathinfo(cigar_str)
            
            # (a) name = input sequence name
            out.write("\t".join([name, paths_str, cigar_str, r_span_str, q_span_str]) + "\n")
    
    
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="program distract overlap cigar and alignments on contigs")
    parser.add_argument("-i", "--input", help="path to input data file",dest="input", type=str, required=True)
    parser.add_argument("-m", "--merge", help="path to input data file",dest="merge", type=str, required=True)
    parser.add_argument("-o", "--output", help="path to input data file",dest="output", type=str, required=True)
    parser.add_argument("-a", "--align", help="path to input data file",dest="align", type=str, required=True )
    parser.set_defaults(func=main)
    args = parser.parse_args()
    args.func(args)
    
