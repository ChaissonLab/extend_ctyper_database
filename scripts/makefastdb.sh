#!/bin/bash

query=$3
database=$2
threads=$1
lcasemask=$4

declare -A seq_lens  # associative array (key: sequence name, value: length)

while IFS= read -r line; do
    if [[ $line == ">"* ]]; then
        # Save previous entry if present
        if [[ -n $name ]]; then
            seq_lens["$name"]=$seqlen
        fi
        # Extract first token from header (before space or tab)
        name=$(echo "${line#>}" | cut -f1 -d' ' | cut -f1)
        seqlen=0
    else
        (( seqlen += ${#line} ))
    fi
done < "$query"

# Save the last entry
if [[ -n $name ]]; then
    seq_lens["$name"]=$seqlen
fi

# Output to file
outfile="${database}.ln"
> "$outfile"  # clear or create file

for name in "${!seq_lens[@]}"; do
    echo -e "$name\t${seq_lens[$name]}" >> "$outfile"
done

lastdb $lcasemask -P $threads $database $query
