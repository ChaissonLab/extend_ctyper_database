#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 input.txt"
  exit 1
fi

infile="$1"
outfile="${infile%.txt}.bed"
mergefile="${infile%.txt}.merge.bed"

awk -F"\t" -v OFS="\t" '
  {
    if (NF < 2) next

    name=$1
    sub(/^>/, "", name)   # strip leading ">"

    pos=$2
    # split on the last colon → left part = contig, right part = coords
    lastc = match(pos, /:[^:]*$/)
    if (lastc == 0) next
    contig = substr(pos, 1, lastc-1)
    coords = substr(pos, lastc+1)

    # split coords on "-" → start, end
    split(coords, b, "-")
    if (length(b) != 2) next
    start=b[1]; end=b[2]

    # strip any trailing + or - from end coord
    sub(/[+-]$/, "", end)

    if (start ~ /^[0-9]+$/ && end ~ /^[0-9]+$/) {
      if (end < start) { t=start; start=end; end=t }
      print contig, start, end, name
    }
  }
' "$infile" > "$outfile"

sort -k1,1 -k2,2n "$outfile" | bedtools merge -i - -c 4 -o distinct > "$mergefile"

totalsize=$(awk '{sum += $3 - $2} END {print sum}' "$mergefile")

# Print filename and total size
echo -e "$(basename "$infile")\t$totalsize"
