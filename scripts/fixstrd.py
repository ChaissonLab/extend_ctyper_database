#!/usr/bin/env python3
"""
Strand-fix pipeline

Usage:
  fixstrd.py -i /path/to/input.fa

Reads:
  1) input.fa
  2) input.fa_graph.FA
  3) input.fa_allgraphalign.out

Writes:
  1) input.fa_graph.FA_strdfix.fa
  2) input.fa_allgraphalign.out_strdfix.fa

Semantics:
  - Input FASTA: name = first field in header (before any whitespace / tab), without '>'.

  - Graph FASTA headers: e.g. CHM13_h1_9_131488_131941
      seqname = first 3 '_' fields (CHM13_h1_9)
      start,end = last 2 fields (0-based, right-open)
      Compare seq with:
        fseq = inputseqs[seqname][start:end]
        rseq = reverse_complement(inputseqs[seqname])[start:end]
      If matches fseq -> '+'
      If matches rseq -> '-'
      If '-', change sequence to reverse complement and header to:
         seqname_{len(inputseqs[seqname]) - end}_{len(inputseqs[seqname]) - start}
      Record FA_header_changed[old] = new

  - Align file (input_allgraphalign.out):
      col1: name in input.fa (same as seqname)
      col2: segments, e.g. "<H1>H2<..."
      col3: segments + CIGAR, e.g. "<H1:CIGAR>H2:CIGAR..."
      col5: query coordinates for each segment, "start_end;start_end;..."

    Step 1 (strand detection):
      - Using ORIGINAL FA_seq (from input_graph.FA) and ORIGINAL col3:
        decode oldseq from col3 CIGAR + FA_seq (full length).
        Insertions and substitutions are represented as 'N's in oldseq.
        Compare oldseq to:
          inputseqs[seqname][:k]                       (forward)
          reverse_complement(inputseqs[seqname])[:k]   (reverse)
        with 'N' matching any base.

    Step 2 (segment/header update):
      - Update names in col2/col3 using FA_header_changed.
      - For each segment whose header changed, flip its local strand char '<'/'>'.

    Step 3 (if line is reverse vs input):
      - If forward: done.
      - If reverse:
          1) Reverse segment order and flip '<'/'>' for all segments in col2 and col3.
          2) Reverse each segment CIGAR and reverse-complement insertion/substitution bases.
          3) Reverse order of col5 segments, and for each "start_end"
             convert to "{len(oldseq) - end}_{len(oldseq) - start}".
             Here len(oldseq) = len(inputseqs[seqname]).

  - CIGAR unit grammar:
      <unit> = <count><op><seq_optional>
      where:
        count: \d+
        op   : one of = M X D I H (case-insensitive)
        seq  : [A-Za-z]* (may be empty; for I/X/S this is the query bases)
"""

import argparse
import gzip
import os
import sys
from typing import Dict, List, Tuple, Optional

DNA_COMP = str.maketrans(
    "ACGTRYMKBDHVNacgtrymkbdhvn",
    "TGCAYRKMVHDBNtgcaYRKMVHDBN"
)


def open_text(path: str, mode: str = "rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


def die(msg: str) -> None:
    print(msg, file=sys.stderr)
    sys.exit(1)


def revcomp(seq: str) -> str:
    return seq.translate(DNA_COMP)[::-1]


def read_fasta_input(path: str) -> Dict[str, str]:
    """
    Input fasta:
      - name = first field in header (split on tab/whitespace), without '>'
      - sequence stored uppercased
    """
    seqs: Dict[str, List[str]] = {}
    name: Optional[str] = None
    with open_text(path, "rt") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                hdr = line[1:].strip()
                first_tab_field = hdr.split("\t", 1)[0]
                name_token = first_tab_field.split()[0]
                if not name_token:
                    raise RuntimeError(f"Empty FASTA header in {path}: {line}")
                name = name_token
                seqs.setdefault(name, [])
            else:
                if name is None:
                    raise RuntimeError(f"Sequence line before header in {path}")
                seqs[name].append(line.strip())
    return {k: "".join(v).upper() for k, v in seqs.items()}


def iter_fasta_records(path: str):
    """
    For graph fasta: header = first whitespace token after '>'.
    Yield (header, seq_upper).
    """
    header = None
    chunks: List[str] = []
    with open_text(path, "rt") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks).upper()
                header = line[1:].split()[0]
                chunks = []
            else:
                if header is None:
                    raise RuntimeError(f"Sequence before header in {path}")
                chunks.append(line.strip())
        if header is not None:
            yield header, "".join(chunks).upper()


def write_fasta(records: List[Tuple[str, str]], path: str, width: int = 60) -> None:
    with open(path, "wt") as w:
        for hdr, seq in records:
            w.write(f">{hdr}\n")
            for i in range(0, len(seq), width):
                w.write(seq[i:i + width] + "\n")


def flip_dir(ch: str) -> str:
    if ch == ">":
        return "<"
    if ch == "<":
        return ">"
    raise ValueError(f"Invalid dir char: {ch}")


def parse_graph_header(h: str) -> Tuple[str, int, int]:
    """
    Graph header like CHM13_h1_9_131488_131941
      seqname = first 3 '_' fields, e.g. CHM13_h1_9
      start,end = last 2 fields
    """
    parts = h.split("_")
    if len(parts) < 5:
        raise RuntimeError(f"Bad graph header (need >=5 '_' fields): {h}")
    seqname = "_".join(parts[:3])
    try:
        start = int(parts[-2])
        end = int(parts[-1])
    except ValueError:
        raise RuntimeError(f"Bad start/end in graph header: {h}")
    return seqname, start, end


def parse_segments_field2(s: str) -> List[Tuple[str, str]]:
    """
    field2: concatenation of (<|>)NAME blocks, e.g. "<A>B<C"
    """
    out: List[Tuple[str, str]] = []
    i = 0
    n = len(s)
    while i < n:
        d = s[i]
        if d not in "<>":
            raise RuntimeError(f"Field2 parse error at pos {i}: {s}")
        j = i + 1
        while j < n and s[j] not in "<>":
            j += 1
        name = s[i + 1:j]
        if not name:
            raise RuntimeError(f"Empty segment name in field2: {s}")
        out.append((d, name))
        i = j
    return out


def parse_segments_field3(s: str) -> List[Tuple[str, str, str]]:
    """
    field3: concatenation of (<|>)NAME:CIGAR blocks; cigar continues until next < or >
    """
    out: List[Tuple[str, str, str]] = []
    i = 0
    n = len(s)
    while i < n:
        d = s[i]
        if d not in "<>":
            raise RuntimeError(f"Field3 parse error at pos {i}: {s}")
        j = i + 1
        while j < n and s[j] != ":":
            j += 1
        if j >= n:
            raise RuntimeError(f"Field3 missing ':' after name: {s}")
        name = s[i + 1:j]
        if not name:
            raise RuntimeError(f"Empty segment name in field3: {s}")
        k = j + 1
        while k < n and s[k] not in "<>":
            k += 1
        cigar = s[j + 1:k].strip()
        if cigar == "":
            raise RuntimeError(f"Empty CIGAR for segment {name}: {s}")
        out.append((d, name, cigar))
        i = k
    return out


def build_field2(segs: List[Tuple[str, str]]) -> str:
    return "".join(d + nm for d, nm in segs)


def build_field3(segs: List[Tuple[str, str, str]]) -> str:
    return "".join(d + nm + ":" + cg for d, nm, cg in segs)


def parse_cigar_ops(cg: str) -> List[Tuple[int, str, Optional[str]]]:
    """
    Parse CIGAR units of the form:
      <count><op><seq_optional>

    where:
      - count is \d+
      - op is one of = M X D I H (case-insensitive)
      - seq_optional is any non-digit letters [A-Za-z]* (may be empty)
    """
    ops: List[Tuple[int, str, Optional[str]]] = []
    i = 0
    n = len(cg)
    while i < n:
        if not cg[i].isdigit():
            raise RuntimeError(f"CIGAR parse error: expected digit at pos {i} in {cg}")
        j = i
        while j < n and cg[j].isdigit():
            j += 1
        count = int(cg[i:j])
        if j >= n:
            raise RuntimeError(f"CIGAR parse error: missing op after count in {cg}")
        op = cg[j]
        j += 1
        k = j
        while k < n and not cg[k].isdigit():
            k += 1
        seq = cg[j:k]
        if seq == "":
            seq = None
        ops.append((count, op, seq))
        i = k
    return ops


def reverse_cigar(cg: str) -> str:
    """
    Reverse operation order. For I/X/S-like (query bases in cigar), reverse-complement the embedded sequence.
    (H/h and D/d have no sequence; just reversed in order like =/M.)
    """
    ops = parse_cigar_ops(cg)
    ops.reverse()
    out: List[str] = []
    for count, op, seq in ops:
        if seq is None:
            out.append(f"{count}{op}")
        else:
            out.append(f"{count}{op}{revcomp(seq.upper())}")
    return "".join(out)


def decode_oldseq_prefix(
    segs3: List[Tuple[str, str, str]],
    segseqs: Dict[str, str],
    want: int = 1000
    ) -> str:
    """
    Reconstruct query(oldseq) prefix using CIGAR + oriented ORIGINAL segment sequence.
    Only needs first `want` bp.

    Semantics:
      '='/'M' : consume ref, emit ref bases
      'X/x'   : consume ref, emit 'N' * count    (ignore explicit seq)
      'I/i'   : emit 'N' * count                 (ignore explicit seq)
      'S/s'   : emit 'N' * count                 (ignore explicit seq)
      'D/d'   : consume ref only
      'H/h'   : hard-clip on reference (consume ref only)
    """
    out: List[str] = []
    out_len = 0

    for d, name, cg in segs3:
        if out_len >= want:
            break
        if name not in segseqs:
            raise RuntimeError(f"Segment sequence not found for {name}")
        base_seq = segseqs[name]
        oriented = base_seq if d == ">" else revcomp(base_seq)
        ref_i = 0

        for count, op, seq in parse_cigar_ops(cg):
            if out_len >= want:
                break

            op_u = op.upper()

            if op_u in ("=", "M"):
                # match: consume reference, emit same bases
                take = oriented[ref_i:ref_i + count]
                if len(take) != count:
                    raise RuntimeError(
                        f"CIGAR '{op}' overruns segment {name}: need {count} at {ref_i}, len={len(oriented)}"
                    )
                out.append(take)
                out_len += count
                ref_i += count

            elif op_u == "X":
                # substitution: consume ref, but emit Ns as wildcard
                out.append("N" * count)
                out_len += count
                ref_i += count

            elif op_u in ("I", "S"):
                # insertion/softclip on query: emit Ns as wildcard, no ref
                out.append("N" * count)
                out_len += count

            elif op_u == "D":
                # deletion: consume ref only
                ref_i += count

            elif op_u == "H":
                # hard clip on reference: consume ref only
                ref_i += count

            else:
                raise RuntimeError(f"Unsupported CIGAR op '{op}' in {cg}")

    return "".join(out)[:want].upper()


def matches_with_N(a: str, b: str) -> bool:
    """
    Compare two sequences of equal length, treating 'N' (in a) as wildcard that matches any base.
    """
    if len(a) != len(b):
        return False
    for ca, cb in zip(a.upper(), b.upper()):
        if ca == 'N':
            continue
        if ca != cb:
            return False
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", required=True, help="Input FASTA path (used as prefix for derived files).")
    args = ap.parse_args()

    in_fa = args.input
    graph_fa = in_fa + "_graph.FA"
    graph_out = in_fa + "_graph.FA_strdfix.fa"
    align_in = in_fa + "_allgraphalign.out"
    align_out = in_fa + "_allgraphalign.out_strdfix.fa"

    if not os.path.exists(in_fa):
        die(f"Error: input fasta not found: {in_fa}")
    if not os.path.exists(graph_fa):
        die(f"Error: graph fasta not found: {graph_fa}")
    if not os.path.exists(align_in):
        die(f"Error: align file not found: {align_in}")

    # Step 1: input fasta into memory
    inputseqs = read_fasta_input(in_fa)

    # Step 2: graph fasta validation and fixing
    header_changed: Dict[str, str] = {}     # old -> new
    segseqs_orig: Dict[str, str] = {}       # ORIGINAL header -> ORIGINAL sequence
    graph_records_out: List[Tuple[str, str]] = []
    seen_headers = set()

    for hdr, seq in iter_fasta_records(graph_fa):
        # store original seq for strand detection
        segseqs_orig[hdr] = seq.upper()

        seqname, start, end = parse_graph_header(hdr)
        if seqname not in inputseqs:
            die(f"Error {hdr}: seqname '{seqname}' not found in input fasta")

        ref = inputseqs[seqname]
        L = len(ref)
        if not (0 <= start <= end <= L):
            die(f"Error {hdr}: coordinates out of range for {seqname} (len={L}): {start}-{end}")

        fseg = ref[start:end].upper()
        rseg = revcomp(ref[L - end:L - start].upper())

        if seq.upper() == fseg:
            out_hdr = hdr
            out_seq = seq.upper()
        elif seq.upper() == rseg:
            new_start = L - end
            new_end = L - start
            out_hdr = f"{seqname}_{new_start}_{new_end}"
            out_seq = revcomp(seq.upper())
            header_changed[hdr] = out_hdr
        else:
            die(f"Error {hdr}: segment does not match forward or reverse-complement slice of {seqname}")

        if out_hdr in seen_headers:
            die(f"Error: duplicate output graph header after fixing: {out_hdr}")
        seen_headers.add(out_hdr)

        graph_records_out.append((out_hdr, out_seq))

    write_fasta(graph_records_out, graph_out)

    # Step 3: align file fixing
    with open(align_in, "rt") as r, open(align_out, "wt") as w:
        for line_no, line in enumerate(r, 1):
            line = line.rstrip("\n")
            if not line:
                w.write("\n")
                continue

            cols = line.split("\t")
            if len(cols) < 5:
                die(f"Error line {line_no}: expected >=5 columns, got {len(cols)}")

            seqname = cols[0]
            if seqname not in inputseqs:
                die(f"Error line {line_no}: col1 '{seqname}' not found in input fasta")

            field2 = cols[1]
            field3 = cols[2]
            field5 = cols[4]

            # --- Step 1: strand detection using ORIGINAL FA_seq and ORIGINAL field3 ---
            segs2_orig = parse_segments_field2(field2)
            segs3_orig = parse_segments_field3(field3)

            ref = inputseqs[seqname]
            try:
                oldseq = decode_oldseq_prefix(segs3_orig, segseqs_orig, want=len(ref))
            except RuntimeError as e:
                print(f"Error format: {align_in}: {seqname}")
                continue

            k = len(oldseq)
            if k == 0:
                die(f"Error line {line_no}: decoded oldseq prefix is empty")

            fseq = ref[:k].upper()
            rseq = revcomp(ref).upper()[:k]

            if matches_with_N(oldseq, fseq):
                line_strand = "+"
            elif matches_with_N(oldseq, rseq):
                line_strand = "-"
            else:
                print(
                    f"Error cigar: {align_in} {line_no}: oldseq does not match forward or reverse of {seqname}",
                )
                continue

            oldseq_len = len(ref)

            # --- Step 2: apply header_changed + per-segment strand flip for changed segments ---
            segs2_u: List[Tuple[str, str]] = []
            for d, nm in segs2_orig:
                if nm in header_changed:
                    segs2_u.append((flip_dir(d), header_changed[nm]))
                else:
                    segs2_u.append((d, nm))
            field2_u = build_field2(segs2_u)

            segs3_u: List[Tuple[str, str, str]] = []
            for d, nm, cg in segs3_orig:
                if nm in header_changed:
                    segs3_u.append((flip_dir(d), header_changed[nm], cg))
                else:
                    segs3_u.append((d, nm, cg))
            field3_u = build_field3(segs3_u)

            # --- Step 3: if line is reverse vs input, flip whole line ---
            if line_strand == "-":
                # reverse segment order, flip dir chars, and reverse CIGARs
                segs2_r = [(flip_dir(d), nm) for (d, nm) in reversed(segs2_u)]
                segs3_r: List[Tuple[str, str, str]] = []
                for d, nm, cg in reversed(segs3_u):
                    segs3_r.append((flip_dir(d), nm, reverse_cigar(cg)))

                field2_u = build_field2(segs2_r)
                field3_u = build_field3(segs3_r)

                # fix col5 (query coordinates per segment on oldseq)
                parts = field5.split(";") if field5 else []
                parts = list(reversed(parts))
                new_parts: List[str] = []
                for p in parts:
                    p = p.strip()
                    if not p:
                        new_parts.append(p)
                        continue
                    if "_" not in p:
                        die(f"Error line {line_no}: bad coord token (expect start_end): '{p}'")
                    a, b = p.split("_", 1)
                    try:
                        st = int(a)
                        en = int(b)
                    except ValueError:
                        die(f"Error line {line_no}: bad coord ints in '{p}'")
                    new_st = oldseq_len - en
                    new_en = oldseq_len - st
                    new_parts.append(f"{new_st}_{new_en}")
                field5_u = ";".join(new_parts)
            else:
                field5_u = field5

            # Write back updated columns (1,2,3,5) preserving others
            cols[1] = field2_u
            cols[2] = field3_u
            cols[4] = field5_u
            w.write("\t".join(cols) + "\n")


if __name__ == "__main__":
    main()

