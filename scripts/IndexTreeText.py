#!/usr/bin/env python3
# IndexTreeText.py

import argparse
from collections import defaultdict
import sys

nametochr = {
    'NC_060925.1': 'chr1',
    'NC_060926.1': 'chr2',
    'NC_060927.1': 'chr3',
    'NC_060928.1': 'chr4',
    'NC_060929.1': 'chr5',
    'NC_060930.1': 'chr6',
    'NC_060931.1': 'chr7',
    'NC_060932.1': 'chr8',
    'NC_060933.1': 'chr9',
    'NC_060934.1': 'chr10',
    'NC_060935.1': 'chr11',
    'NC_060936.1': 'chr12',
    'NC_060937.1': 'chr13',
    'NC_060938.1': 'chr14',
    'NC_060939.1': 'chr15',
    'NC_060940.1': 'chr16',
    'NC_060941.1': 'chr17',
    'NC_060942.1': 'chr18',
    'NC_060943.1': 'chr19',
    'NC_060944.1': 'chr20',
    'NC_060945.1': 'chr21',
    'NC_060946.1': 'chr22',
    'NC_060947.1': 'chrX',
    'NC_060948.1': 'chrY',
    }


def normalize_token(s: str) -> str:
    """
    Normalize a token to match tree leaf parsing:
      - strip whitespace
      - remove quotes
      - take first field (split by whitespace/tab)
      - drop anything after ':' (distance safety)
      - drop anything after ',' (for multi-allele cells)
    """
    if s is None:
        return ""
    s = s.strip().replace("'", "").replace('"', "")
    if not s:
        return ""
    s = s.split()[0].split("\t")[0]
    if ":" in s:
        s = s.split(":", 1)[0]
    if "," in s:
        s = s.split(",", 1)[0]
    return s


class node:
    def __init__(self, parent=None, name="", distance=0.0, index=-1):
        self.children = []
        self.parent = parent
        self.name = name
        self.distance = distance
        self.index = index
        self.index2 = index
        self.annotation = (0, 0)
        if parent is not None:
            parent.children.append(self)

    def __str__(self):
        if len(self.children) == 0:
            return self.name + ":" + "{:.7f}".format(self.distance)
        return "(" + ",".join([str(x) for x in self.children]) + "):" + "{:.7f}".format(self.distance)

    def push(self, name="", distance=0.0, index=0):
        return node(self, name, distance, index)

    def str(self):
        return str(self)

    def build(self, text, allnodes=None):
        """
        Parse a Newick string into this node tree.
        Compatible with DistractTree.py.
        """
        if allnodes is None:
            allnodes = []
        if not text:
            return self
        if text.endswith(";"):
            text = text[:-1]

        allnodes.append(self)
        current_node = self
        index = 0
        allnames = []

        for char in text:
            if char in [" ", "'"]:
                continue

            if char == "(":
                current_node = current_node.push()
                allnodes.append(current_node)

            elif char == ")":
                name, distance = (current_node.name.split(":") + ["0.0"])[:2]
                name = normalize_token(name)
                if name:
                    current_node.name = name
                    allnames.append(name)
                    if len(current_node.children) == 0:
                        current_node.index = index
                        index += 1
                else:
                    allnames.append(current_node.name)

                try:
                    current_node.distance = float(distance.strip())
                except Exception:
                    current_node.distance = 0.0

                current_node = current_node.parent
                if current_node is not None and len(current_node.name) == 0:
                    current_node.name = "bh_" + str(len(allnames))

            elif char == "," or char == ";":
                name, distance = (current_node.name.split(":") + ["0.0"])[:2]
                name = normalize_token(name)
                if name:
                    current_node.name = name
                    allnames.append(name)
                    if len(current_node.children) == 0:
                        current_node.index = index
                        index += 1
                else:
                    allnames.append(current_node.name)

                try:
                    current_node.distance = float(distance.strip())
                except Exception:
                    current_node.distance = 0.0

                if current_node.parent is not None:
                    current_node = current_node.parent.push()
                    allnodes.append(current_node)

            else:
                current_node.name += char

        return self


def build_leaf_index_map(row_line: str) -> dict:
    """
    From a matrix row line, build a mapping:
      leaf_name -> index

    Rule for indices (1-based columns):
      - col1   -> index 1
      - col2-4 -> ignored
      - col5   -> index 2
      - col6   -> index 3
      - ...    -> index = col_num - 3 for col_num >= 5

    We use normalize_token on each cell, and only map non-empty tokens.
    If multiple columns share the same token, the first column wins.
    """
    cols = row_line.split("\t")
    leaf_to_index = {}

    for i, col in enumerate(cols):
        colnum = i + 1  # 1-based column number

        # Skip 2nd, 3rd, 4th columns entirely
        if colnum in (2, 3, 4):
            continue

        token = normalize_token(col)
        if not token:
            continue

        if colnum == 1:
            idx = 1
        else:
            # colnum >= 5
            idx = colnum - 3

        if token not in leaf_to_index:
            leaf_to_index[token] = idx

    return leaf_to_index


def chrom_key(chrom: str):
    """
    Produce a sortable key for chromosome names like 'chr1'..'chr22','chrX','chrY'.
    Numeric chromosomes come first in numeric order; X,Y,M (if present) after.
    Any unknown names are put at the end, lexicographically.
    """
    c = chrom
    if c.startswith("chr"):
        c = c[3:]
    # numeric
    if c.isdigit():
        return (0, int(c))
    # common non-numeric
    if c in ("X", "x"):
        return (1, 23)
    if c in ("Y", "y"):
        return (1, 24)
    if c in ("M", "MT", "m", "mt"):
        return (1, 25)
    # fallback: treat as string
    return (2, c)


def main():
    ap = argparse.ArgumentParser(
        description="Replace tree leaf names with column indices, map CHM13 contigs to chr names, and sort by chrom,start."
    )
    ap.add_argument("-i", "--input", required=True,
                    help="Input file: output of DistractTree.py (with #$ tree lines and matrix lines).")
    ap.add_argument("-o", "--output", required=True,
                    help="Output file with indexed tree leaf names, remapped chrom column, and sorted rows.")
    ap.add_argument("--debug", action="store_true")
    args = ap.parse_args()

    header_lines = []
    records = []  # list of (sort_key, tree_line_str, row_line_str)

    with open(args.input, "rt") as fin:
        # 1) Read header block (everything before first '#$' line)
        while True:
            pos = fin.tell()
            line = fin.readline()
            if not line:
                break
            if line.startswith("#$"):
                # rewind to start of this tree line
                fin.seek(pos)
                break
            header_lines.append(line.rstrip("\n"))

        # 2) Read tree+row pairs
        while True:
            tree_line = fin.readline()
            if not tree_line:
                break  # EOF

            tree_line = tree_line.rstrip("\n")
            if not tree_line.startswith("#$"):
                # Unexpected extra non-tree line after header; keep as extra header/footer
                header_lines.append(tree_line)
                continue

            row_line = fin.readline()
            if not row_line:
                # No accompanying row line; store tree only as-is (no sorting key) at end
                header_lines.append(tree_line)
                break
            row_line = row_line.rstrip("\n")

            # Split row into columns (tab-separated), change 2nd column via nametochr
            cols = row_line.split("\t")
            if len(cols) >= 2:
                orig_contig = cols[1]
                cols[1] = nametochr.get(orig_contig, orig_contig)
            # Rebuild row_line with updated chrom name in col2
            row_line = "\t".join(cols)

            # Build mapping from leaf name to index using the updated matrix row
            leaf_to_index = build_leaf_index_map(row_line)

            # Parse tree
            newick_text = tree_line[2:]  # strip leading "#$"
            newick_text = newick_text.strip()
            if newick_text.endswith(";"):
                newick_text = newick_text[:-1]

            root = node()
            allnodes = []
            root.build(newick_text, allnodes=allnodes)

            # Replace leaf names with indices when possible
            for nd in allnodes:
                if len(nd.children) == 0 and nd.name:
                    tok = normalize_token(nd.name)
                    idx = leaf_to_index.get(tok)
                    if idx is not None:
                        nd.name = str(idx)
                    else:
                        if args.debug:
                            sys.stderr.write(f"[warn] leaf '{nd.name}' (norm '{tok}') has no column index\n")

            new_tree_line = "#$" + root.str() + ";"

            # Build sorting key from 2nd and 3rd columns (chrom, start)
            chrom = cols[1] if len(cols) > 1 else ""
            try:
                start = int(cols[2]) if len(cols) > 2 else 0
            except ValueError:
                start = 0

            sort_key = (chrom_key(chrom), start)

            records.append((sort_key, new_tree_line, row_line))

    # 3) Sort records by chrom,start
    records.sort(key=lambda x: x[0])

    # 4) Write output: header (unsorted) then sorted records
    with open(args.output, "wt") as out:
        for h in header_lines:
            out.write(h + "\n")
        for _, tree_line_str, row_line_str in records:
            out.write(tree_line_str + "\n")
            out.write(row_line_str + "\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

