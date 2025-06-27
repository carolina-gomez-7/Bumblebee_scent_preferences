#!/bin/bash
##################################################################################
# Script:     04_cleanHeaders.sh
# Author:     Carolina Gomez Ramirez & Maximiliano Zuluaga Forero
# Purpose:    Clean header names from gb or fasta files keeping only the organism field
# Usage: ./04_cleanHeaders.sh input_file.fasta
##################################################################################


if [ -z "$1" ]; then
  echo "Usage: $0 <input_fasta_file>"
  exit 1
fi

input="$1"
base=$(basename "$input")
output="${base%.*}_clean.${base##*.}"

gawk '/^>/ {
  match($0, /\[organism=([^\]]+)\]/, arr);
  if (arr[1]) {
    gsub(" ", "_", arr[1]);
    print ">" arr[1];
  } else {
    print $0
  }
  next
} { print }' "$input" > "$output"

echo "Cleaned file saved as: $output"


