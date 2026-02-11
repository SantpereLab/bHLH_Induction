#!/bin/bash

cd ~/proneural/chip-seq
# Define the list of accessions
# Read the list of accessions from the file
accessions=( $(cat srx_accessions.txt) )


# Iterate over each accession and download the peaks. For human
for accession in "${accessions[@]}"; do
    echo $accession
    url="https://chip-atlas.dbcls.jp/data/hg38/eachData/bed05/${accession}.05.bed"
    sbatch --wrap="wget -O \"${accession}.bed\" \"$url\""
done


# Iterate over each accession and download the peaks for mm10 if the file is empty
for accession in "${accessions[@]}"; do
    if [ ! -s "${accession}.bed" ]; then
        echo "Downloading ${accession} for mm10"
        url="https://chip-atlas.dbcls.jp/data/mm10/eachData/bed05/${accession}.05.bed"
        sbatch --wrap="wget -O \"${accession}.bed\" \"$url\""
    fi
done


# Create an associative array to map SRX to GSM
declare -A srx_to_gsm

# Read the mapping file and populate the associative array
while read -r gsm srx; do
    srx_to_gsm[$srx]=$gsm
done < srx_accessions_df.txt

# Rename the .bed files from SRX accession to GSM accession
for accession in "${accessions[@]}"; do
    if [ -f "${accession}.bed" ]; then
        gsm_accession=${srx_to_gsm[$accession]}
        if [ -n "$gsm_accession" ]; then
            mv "${accession}.bed" "${gsm_accession}.bed"
        else
            echo "No GSM accession found for ${accession}"
        fi
    fi
done


# Check which accessions have not been successfully downloaded (empty files)
for accession in "${accessions[@]}"; do
    if [ ! -s "${srx_to_gsm[$accession]}.bed" ]; then
        echo "Failed to download or empty file: ${accession}"
    fi
done

