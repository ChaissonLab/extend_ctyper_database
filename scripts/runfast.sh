#!/bin/bash

query=$1
database=$2
threads=$3
simi="${4:-90}"
size="${5:-300}"

dblenfile=$database".ln"

lastal  -j 3 -e 300 -d 50  -P $threads $database $query  | maf-convert sam | \
awk -v cutoff=$size -v simicutoff=$simi -v lenfile=$dblenfile '
BEGIN {
    OFS="\t"
    # Load the length file ONCE
    while ((getline < lenfile) > 0) {
        seq_lens[$1] = $2
    }
    close(lenfile)
} 
$1 ~ /^@/ {print; next} 
{
            
    # Determine strand: if leading H, assume reverse
    strand = ($2 == 16) ? -1 : 1
    cigar = $6

    match(cigar, /^([0-9]+)H/, a); hq_1 = (a[1] ? a[1] : 0)
    match(cigar, /([0-9]+)H$/, b); hq_2 = (b[1] ? b[1] : 0)

    newstart = hq_1 + 1
    if (strand == -1){
        newstart = hq_2 + 1
    } 

    sub(/^[0-9]+H/, "", cigar)
    sub(/[0-9]+H$/, "", cigar)

    gsub(/=/, "M", cigar)
    gsub(/I/, "_", cigar)
    gsub(/D/, "I", cigar)
    gsub(/_/, "D", cigar)
    cindex=0 
    len_aln=0 
    len_q = 0 
    len_r = 0
    oldcigar = cigar
    while (match(cigar, /[0-9]+[M=XIDH]/)) {
        n = substr(cigar, RSTART, RLENGTH - 1)  # Number
        op = substr(cigar, RSTART + length(n), 1)  # Operation
    
        # Update alignment length
        if (op ~ /[M=XID]/) len_aln += n
        if (op ~ /[M=XI]/)  len_r += n   
        if (op ~ /[M=XD]/) len_q += n    

        if (strand == -1) {
            cindex++
            tok[cindex] = n op
        }
    
        # Remove parsed part from cigar
        cigar = substr(cigar, RSTART + RLENGTH)
    }
    
    
    if (strand == -1) {
        cigar = ""
        for (i = cindex; i >= 1; i--) {
            cigar = cigar tok[i]
        }
    } else {
        cigar = oldcigar
    }

    hr_1 = $4 - 1
    hr_2 = seq_lens[$3] - $4 - len_r + 1  # safer calculation

    hr_1_s = (hr_1 > 0) ? hr_1 "H" : ""
    hr_2_s = (hr_2 > 0) ? hr_2 "H" : ""

    if (strand == -1) {
        cigar = hr_2_s cigar hr_1_s
    } else {
        cigar = hr_1_s cigar hr_2_s
    }


    # Parse NM:i: (mismatches/indels)
    split($12, nm_tag, ":")
    nm = nm_tag[3]
    if (len_aln > cutoff ) {
        # Calculate identity score (based on mismatches/indels and alignment length)
        idscore = (len_aln > 0) ? (100*(1 - nm / len_aln)) : 0
        # If identity score passes threshold
        if (idscore > simicutoff) {
            # Replace SEQ field with "*"
            $10 = "*"
            $4 = newstart
            tmp = $1; $1 = $3; $3 = tmp
            # Replace NM:i:, AS:i:, and append PI:f:
            $12 = "NM:i:" nm
            $13 = "AS:i:" len_aln  # Use the alignment length for AS (as a proxy)
            $15 = "PI:f:" idscore  # Set the computed PI
            $6 = cigar
            print 
            # Print the updated SAM line
        }
    }
}'

