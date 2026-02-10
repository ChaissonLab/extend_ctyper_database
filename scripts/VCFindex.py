#!/usr/bin/env python3

import argparse
import collections as cl
import pysam

BIN = 10000  # coordinate bin size for indexing regions (tune if needed)


def to_ranges(nums):
    if not nums:
        return []
    nums = sorted(nums)
    ranges = []
    start = end = nums[0]
    for n in nums[1:]:
        if n == end + 1:
            end = n
        elif n > end + 1:
            ranges.append((start, end))
            start = end = n
    ranges.append((start, end))
    return ranges


def sample_key_from_name(name: str) -> str:
    parts = name.split("_")
    if len(parts) >= 3:
        return "_".join(parts[1:3])
    return ""


def prefix_from_name(name: str) -> str:
    # "text before '_'"
    p = name.split("_", 1)
    return p[0] if p else ""


def interval_dist(a0, a1, b0, b1) -> int:
    """Distance between half-open intervals [a0,a1) and [b0,b1). Overlap -> 0"""
    if a1 <= b0:
        return b0 - a1
    if b1 <= a0:
        return a0 - b1
    return 0


def interval_olen(a0, a1, b0, b1) -> int:
    """Overlap length of half-open intervals."""
    x0 = a0 if a0 > b0 else b0
    x1 = a1 if a1 < b1 else b1
    return x1 - x0 if x1 > x0 else 0


def parse_info_end_1based(info: str):
    # END=#### in VCF is typically 1-based inclusive.
    # We'll convert to 0-based exclusive by using end0 = END (works for POS=END single-base case).
    for tok in info.split(";"):
        if tok.startswith("END="):
            v = tok[4:]
            try:
                return int(v)
            except ValueError:
                return None
    return None


def load_regions_and_refprefixes(qregionfile: str, samplename: str, refkey: str = "HG38_h1"):
    """
    region TSV columns:
        name chrom strand start end thickstart thickend  (your earlier)
    BUT your current code expects:
        name  contig  strd  start  end  left  right
    We'll keep your existing semantics: start/end + left/right => sub_start/sub_end

    This loader returns:
        regions_by_contig: only records matching samplename
        refprefix_by_contig: reference intervals for each prefix from records with sample_key == refkey
                            { contig: { prefix: [(s,e), ...] } }
    """
    regions_by_contig = cl.defaultdict(list)
    refprefix_by_contig = cl.defaultdict(lambda: cl.defaultdict(list))
    
    with open(qregionfile, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            # must have at least 7 cols because we read parts[6]
            if len(parts) < 8:
                continue
            
            contig = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            name = parts[3].lstrip(">")
            score = int(parts[4])
            strd = parts[5]
            left = int(parts[6])
            right = int(parts[7])
            
            skey = sample_key_from_name(name)
            pref = prefix_from_name(name)
            
            valid = (left != -1 and right != -1)
            sub_start = start + left if valid else -1
            sub_end = end - right if valid else -1
            if valid and sub_start >= sub_end:
                valid = False
                sub_start = sub_end = -1
                
            rec = {
                "name": name,
                "prefix": pref,
                "contig": contig,
                "score": score,
                "strd": strd,
                "start": start,
                "end": end,
                "left": left,
                "right": right,
                "valid": valid,
                "sub_start": sub_start,
                "sub_end": sub_end,
            }
            
            # A) record HG38 coords per prefix from HG38_h1 rows
            if skey == refkey:
                if valid:
                    ref_s, ref_e = sub_start, sub_end
                else:
                    ref_s, ref_e = start, end
                if ref_s < ref_e:
                    refprefix_by_contig[contig][pref].append((ref_s, ref_e))
                    
            # sample-specific query regions
            if skey == samplename:
                regions_by_contig[contig].append(rec)
                
    return regions_by_contig, refprefix_by_contig


def build_region_bin_index(records, slop: int):
    """
    Bin-map for valid subregions (expanded by slop) on one contig:
        bin_id -> [record_index, ...]
    """
    bmap = cl.defaultdict(list)
    for i, rec in enumerate(records):
        if not rec["valid"]:
            continue
        s = rec["sub_start"]
        e = rec["sub_end"]
        if s >= e:
            continue
        s2 = max(0, s - slop)
        e2 = e + slop
        b0 = s2 // BIN
        b1 = (e2 - 1) // BIN if e2 > 0 else b0
        for b in range(b0, b1 + 1):
            bmap[b].append(i)
    return bmap


def choose_prefix_by_ref(refprefix_by_contig, tie_prefixes, ref_contig, r0, r1):
    """
    C/D) Use variant reference interval [r0,r1) to choose prefix:
        - choose prefix with largest overlap size with its HG38 interval(s)
        - if all overlaps are 0, choose prefix with smallest distance to HG38 interval(s)
        - tie -> first in tie_prefixes order (stable)
    """
    prefmap = refprefix_by_contig.get(ref_contig)
    if not prefmap:
        return tie_prefixes[0]
    
    best_pref = tie_prefixes[0]
    best_olen = -1
    best_dist = 10**18
    
    for pref in tie_prefixes:
        ivs = prefmap.get(pref)
        if not ivs:
            continue
        
        # summarize across potentially multiple HG38 intervals for this prefix
        max_olen = 0
        min_dist = 10**18
        for s, e in ivs:
            if s >= e:
                continue
            olen = interval_olen(r0, r1, s, e)
            if olen > max_olen:
                max_olen = olen
            d = interval_dist(r0, r1, s, e)
            if d < min_dist:
                min_dist = d
                
        # prefer larger overlap; then smaller distance; then earlier order
        if (max_olen > best_olen) or (max_olen == best_olen and min_dist < best_dist):
            best_olen = max_olen
            best_dist = min_dist
            best_pref = pref
            
    return best_pref

def compress_assigned_to_byte_ranges(varts, assigned_vidxs):
    """
    varts entries contain:
        [qry_start, qry_end, lindex, voff, vend, ref_contig, ref0, ref1, varid]
    assigned_vidxs: list of indices into varts
    Returns list of (start_voff, end_vend) ranges (virtual offsets).
    """
    if not assigned_vidxs:
        return []
    
    line_indices = sorted(varts[vidx][2] for vidx in assigned_vidxs)
    ranges = to_ranges(line_indices)
    
    lindex_to_vidx = {varts[vidx][2]: vidx for vidx in assigned_vidxs}
    
    out = []
    for r_start, r_end in ranges:
        v_start = lindex_to_vidx[r_start]
        v_end = lindex_to_vidx[r_end]
        out.append((varts[v_start][3], varts[v_end][4]))  # voff..vend
    return out

def vcfindex(inputfile, qregionfile, rregionfile, output, samplename, slop=3000, refkey="HG38_h1"):
    regions_by_contig, refprefix_by_contig = load_regions_and_refprefixes(
        qregionfile, samplename, refkey=refkey
    )
    if not regions_by_contig:
        print(f"Error: sample {samplename} not found in region TSV (or no parsable rows).")
        return
    
    contigs_wanted = set(regions_by_contig.keys())
    vcfcoordi = cl.defaultdict(list)
    
    # Scan VCF once; store query intervals + virtual offsets + reference interval + varid
    with pysam.BGZFile(inputfile, mode="rb") as f:
        lindex = 0
        while True:
            voff = f.tell()
            line = f.readline()
            if not line:
                break
            vend = f.tell()
            
            line_str = line.decode(errors="replace")
            if line_str.startswith("#"):
                continue
            
            cols = line_str.split("\t")
            if len(cols) < 8:
                lindex += 1
                continue
            
            # reference coordinates from VCF columns
            ref_contig = cols[0]
            try:
                pos1 = int(cols[1])  # 1-based
            except ValueError:
                lindex += 1
                continue
            ref0 = pos1 - 1
            info_field = cols[7]
            end1 = parse_info_end_1based(info_field)
            ref1 = end1 if end1 is not None else pos1  # 0-based exclusive
            
            # varid: prefer VCF ID, else CHROM:POS:REF:ALT
            vcf_id = cols[2] if len(cols) > 2 else ""
            if vcf_id and vcf_id != ".":
                varid = vcf_id
            else:
                ref = cols[3] if len(cols) > 3 else ""
                alt = cols[4] if len(cols) > 4 else ""
                varid = f"{ref_contig}:{pos1}:{ref}:{alt}"
                
            # query coordinates from INFO QRY_REGION=ctg:s-e
            region_strs = [x for x in info_field.split(";") if x.startswith("QRY_REGION=")]
            if not region_strs:
                lindex += 1
                continue
            
            region = region_strs[0][len("QRY_REGION="):]
            if ":" not in region or "-" not in region:
                lindex += 1
                continue
            
            qry_contig, coord_str = region.split(":", 1)
            if qry_contig not in contigs_wanted:
                lindex += 1
                continue
            
            try:
                s_str, e_str = coord_str.split("-", 1)
                q0 = int(s_str)
                q1 = int(e_str)
            except ValueError:
                lindex += 1
                continue
            
            vcfcoordi[qry_contig].append([q0, q1, lindex, voff, vend, ref_contig, ref0, ref1, varid])
            lindex += 1
            
    varassign_path = output + ".varassign"
    
    with open(output, "w") as w, open(varassign_path, "w") as wva:
        for qry_contig, records in regions_by_contig.items():
            varts = vcfcoordi.get(qry_contig, [])
            assigned = [[] for _ in records]
            
            assigned_varid_to_prefix = {}
            assigned_varid_order = []
            
            bmap = build_region_bin_index(records, slop)
            
            for vidx, v in enumerate(varts):
                q0, q1 = v[0], v[1]
                varid = v[8]
                
                # collect candidate records within slop (using bin acceleration)
                s2 = max(0, q0 - slop)
                e2 = q1 + slop
                b0 = s2 // BIN
                b1 = (e2 - 1) // BIN if e2 > 0 else b0
                
                cand = []
                seen = set()
                for b in range(b0, b1 + 1):
                    for ridx in bmap.get(b, []):
                        if ridx in seen:
                            continue
                        seen.add(ridx)
                        rec = records[ridx]
                        rs, re = rec["sub_start"], rec["sub_end"]
                        d = interval_dist(q0, q1, rs, re)
                        if d <= slop:
                            cand.append((d, ridx))
                            
                if not cand:
                    continue
                
                # ---- NEW LOGIC ----
                # 1) Choose prefix by reference FIRST (among all feasible prefixes)
                cand_by_prefix = cl.defaultdict(list)  # prefix -> list of (d, ridx)
                for d, ridx in cand:
                    p = records[ridx]["prefix"]
                    if p:
                        cand_by_prefix[p].append((d, ridx))
                        
                if not cand_by_prefix:
                    continue
                
                prefixes = list(cand_by_prefix.keys())
                if len(prefixes) == 1:
                    var_prefix = prefixes[0]
                else:
                    ref_contig, r0, r1 = v[5], v[6], v[7]
                    var_prefix = choose_prefix_by_ref(refprefix_by_contig, prefixes, ref_contig, r0, r1)
                    
                # record varid -> prefix
                if varid not in assigned_varid_to_prefix:
                    assigned_varid_order.append(varid)
                assigned_varid_to_prefix[varid] = var_prefix
                
                # 2) Within that prefix, choose record by min query distance
                lst = cand_by_prefix.get(var_prefix)
                if not lst:
                    continue
                lst.sort(key=lambda x: (x[0], x[1]))  # min d, then smallest ridx
                winner = lst[0][1]
                assigned[winner].append(vidx)
                
            # write per-contig var assignment lines
            for varid in assigned_varid_order:
                wva.write(f"{varid}:{assigned_varid_to_prefix[varid]}\n")
                
            # write indexed ranges (virtual offsets) per record
            for ridx, rec in enumerate(records):
                if not rec["valid"]:
                    result = ""
                else:
                    byte_ranges = compress_assigned_to_byte_ranges(varts, assigned[ridx])
                    result = ";".join([f"{a}_{b}" for a, b in byte_ranges]) if byte_ranges else ""
                    
                w.write(
                    f"{rec['contig']}\t{rec['start']}\t{rec['end']}\t{rec['name']}\t"
                    f"{rec['score']}\t{rec['strd']}\t{rec['left']}\t{rec['right']}\t{result}\n"
                )
                
def main(args):
    vcfindex(args.input, args.qregion, args.rregion,args.output, args.sample, slop=args.slop, refkey=args.refkey)
    
    
def run():
    parser = argparse.ArgumentParser(description="Index VCF BGZF byte ranges per (sub)region with slop + exclusive prefix assignment.")
    parser.add_argument("-i", "--input",  dest="input",  type=str, required=True, help="input .vcf.gz (BGZF/bgzip)")
    parser.add_argument("-r", "--rregion", dest="rregion", type=str, default = "", help="region TSV from liftover")
    parser.add_argument("-q", "--qregion", dest="qregion", type=str, required=True, help="region TSV from first script")
    parser.add_argument("-o", "--output", dest="output", type=str, required=True, help="output TSV")
    parser.add_argument("-s", "--sample", dest="sample", type=str, required=True, help="sample name key (parts[1:3])")
    parser.add_argument("--slop", dest="slop", type=int, default=300, help="allow gap distance (bp) between variant and region (default 300)")
    parser.add_argument("--refkey", dest="refkey", type=str, default="HG38_h1", help="2nd_3rd key used as reference rows (default HG38_h1)")
    parser.set_defaults(func=main)
    args = parser.parse_args()
    args.func(args)
    
    
if __name__ == "__main__":
    run()
