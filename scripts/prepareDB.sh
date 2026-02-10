#!/bin/bash

thread=64

inputblock=$1

scriptfolder="/project2/mchaisso_100/walfred/projects/newphase/scripts/"

kmermaskCHM13="/project/mchaisso_100/cmb-16/walfred/database/CHM13v2.0_kmermask.fa"
winnowCHM13="/project/mchaisso_100/cmb-16/walfred/database/GCF_009914755.1_T2T-CHM13v2.0_genomic.fa_wmask.fa"
genecode="/project/mchaisso_100/cmb-16/walfred/database/catLiftOffGenesV1.gff3"
destinyfolder="/project2/mchaisso_100/walfred/projects/newphase/Build/firstrun/groups/"

mkdir -p blocks groups blocks_exons groups_exons $destinyfolde

# Step 1: kmermask
python "$scriptfolder/getblocks.py" -i "$inputblock" -r "$kmermaskCHM13" -o blocks/
bash "$scriptfolder/makegroup.sh" "$inputblock" blocks/ groups/

ls groups/ -1 | xargs -P $thread -I {} bash -c '
    name="${1%.fa}"
    mkdir -p "$2/$name"
    mv "groups/$name.fa" "$2/$name/"
' _ {} "$destinyfolder"

rm -f blocks/*
rm -r groups/*
# Step 2: winnowmask
python "$scriptfolder/getblocks.py" -i "$inputblock" -r "$winnowCHM13" -o blocks/
bash "$scriptfolder/makegroup.sh" "$inputblock" blocks/ groups/

ls groups/ -1 | xargs -P $thread -I {} bash -c '
    name="${1%.fa}"
    mv "groups/$name.fa" "$2/$name/${name}_origin.fa"
' _ {} "$destinyfolder"

# Step 3: exons
python "$scriptfolder/block_exons.py" -i "$inputblock" -r "$kmermaskCHM13" -g "$genecode" -o blocks_exons/
bash "$scriptfolder/makegroup_exon.sh" "$inputblock" blocks_exons/ groups_exons/

ls groups_exons/ -1 | xargs -P $thread -I {} bash -c '
    name="${1%_exons.fa}"
    mv "groups_exons/${name}_exons.fa" "$2/$name/"
' _ {} "$destinyfolder"

wait


# Step 4: makeblastdb for *_origin.fa
cat $inputblock | cut -f 8 | sort | uniq  | xargs -P $thread -I {} bash -c '
    f="$1/{}"
    if [ -f "$f/{}_origin.fa" ]; then
        makeblastdb -in "$f/{}_origin.fa" -dbtype nucl -parse_seqids -out "$f/{}.fa_db" > /dev/null 2>&1
    fi
' _ "$destinyfolder"

# Step 5: makeblastdb for *_exons.fa
cat $inputblock | cut -f 8 | sort | uniq | xargs -P $thread -I {} bash -c '
    f="$1/{}"
    if [ -f "$f/{}_exons.fa" ]; then
        makeblastdb -in "$f/{}_exons.fa" -dbtype nucl -parse_seqids -out "$f/{}_exons.fa_db" > /dev/null 2>&1
    fi
' _ "$destinyfolder"
