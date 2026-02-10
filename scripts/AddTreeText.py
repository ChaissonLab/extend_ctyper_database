stractTree.py

import argparse
import collections as cl
import sys


def normalize_token(s: str) -> str:
    """
    Normalize a token to match tree leaf parsing:
      - strip whitespace
      - remove quotes
      - take first field (split by whitespace/tab)
      - drop anything after ':' (distance safety)
      - drop anything after ',' (multi-allele cells)
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


def prefix_from_name(name: str) -> str:
    name = normalize_token(name)
    return name.split("_", 1)[0] if "_" in name else name


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

    def filter_by_leafnames_keep(self, keep_names: set):
        # leaf
        if len(self.children) == 0:
            return self if normalize_token(self.name) in keep_names else None

        # internal
        new_children = []
        for ch in self.children:
            kept = ch.filter_by_leafnames_keep(keep_names)
            if kept is not None:
                new_children.append(kept)

        self.children = new_children
        for ch in self.children:
            ch.parent = self

        if len(self.children) == 0:
            return None

        if len(self.children) == 1:
            only = self.children[0]
            only.distance += self.distance
            only.parent = self.parent
            return only

        return self


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--Input", required=True)
    ap.add_argument("-d", "--database", required=True)
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("--debug", action="store_true")
    args = ap.parse_args()

    # prefix -> list of per-line keep sets (each set contains normalized tokens from that line)
    used_lines = cl.defaultdict(list)

    header = ""
    # Read input: for each line, store set(tokens)
    with open(args.Input, "rt") as f:
        index = 0
        for raw in f:
            line = raw.rstrip("\n")
            if index == 0 and line.startswith("#"):
                header = line
                index += 1
                continue

            index += 1
            if not line or line.startswith("#"):
                continue

            cols = line.split()  # tabs or spaces
            if not cols:
                continue

            pref = prefix_from_name(cols[0])

            keep_set = set()
            for tok in cols:
                t = normalize_token(tok)
                if t:
                    keep_set.add(t)

            used_lines[pref].append((keep_set, line))  # keep_set for filtering, and raw line for output

    with open(args.output, "wt") as out, open(args.database, "rt") as fdb:

        # write header (single line), unchanged
        if header:
            out.write(header + "\n")

        for raw in fdb:
            text = raw.strip()
            if not text:
                continue

            # Build tree once for this database line
            root0 = node()
            allnodes0 = []
            root0.build(text, allnodes=allnodes0)

            leaves0 = [n for n in allnodes0 if len(n.children) == 0 and n.name]
            if not leaves0:
                continue

            tree_pref = prefix_from_name(leaves0[0].name)
            cases = used_lines.get(tree_pref)
            if not cases:
                continue

            tree_leaf_set0 = None
            if args.debug:
                tree_leaf_set0 = {normalize_token(x.name) for x in leaves0}
                sys.stderr.write(f"[debug] prefix={tree_pref} cases={len(cases)} tree_leaves={len(tree_leaf_set0)}\n")

            # For each input-line keep set, rebuild+filter (so one case can't mutate another)
            for keep_set, raw_line in cases:
                # rebuild fresh tree (filter mutates nodes)
                root = node()
                root.build(text)

                if args.debug and tree_leaf_set0 is not None:
                    matched = len(keep_set & tree_leaf_set0)
                    sys.stderr.write(f"[debug]   case_keep={len(keep_set)} matched={matched}\n")

                newroot = root.filter_by_leafnames_keep(keep_set)
                if newroot is None:
                    continue
                newroot.parent = None

                # 1) subtree
                out.write("#$" + newroot.str() + ";\n")
                # 2) original input line
                out.write(raw_line + "\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

