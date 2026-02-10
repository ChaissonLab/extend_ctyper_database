#!/bin/bash

if [[ $# -ne 3 ]]; then
    echo "Usage: $0 table_file input_dir output_dir"
    exit 1
fi

table_file=$1
input_dir=$2
output_dir=$3

while read -r line; do
    x=$(echo "$line" | awk '{print $(NF-1)}')  # 2nd last column
    y=$(echo "$line" | awk '{print $NF}')      # Last column

    if [[ -f $input_dir/"$x"_exons.fa ]]; then
        cat $input_dir/"$x"_exons.fa >> $output_dir/"$y"_exons.fa
    else
        echo "Warning: File $input_dir/"$x"_exons.fa not found."
    fi
done < "$table_file"
