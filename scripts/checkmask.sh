#!/usr/bin/env bash

fasta_file="$1"

# Redirect stdout and stderr to temporary named pipes to capture both counts
exec 3> >(read total; echo -n "$total ")
grep -v '^>' "$fasta_file" | tr -d '\n' | tee >(wc -c >&3) | tr -cd 'a-z' | wc -c
