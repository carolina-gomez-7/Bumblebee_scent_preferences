#!/bin/bash

################################################################################
# Script:     01_cleanDatabase_V02.sh
# Author:     Carolina Gomez Ramirez & Maximiliano Zuluaga Forero
# Purpose:    Filter, annotate, and prepare bumblebee Sanger sequence data
# Usage:      ./01_cleanDatabase.sh inputfile.csv
# Output:     inputfile_filtered.csv
#
# REQUIREMENTS:
# 1) The input must be a .csv with the following column structure:
#    ID, Species, Species ID Morphology, Primer, Sequencing code,
#    Nanodrop concentration, Sex, Sequences
#
# 2) Primer names must be normalized to exactly one of:
#    - "Hymeno 1", "Hymeno 2", "Lep F", "Lep R", "16 S IR", "16 S WB"
#
# PROCESSING STEPS:
# 1. Removes carriage returns (\r) for clean line endings.
# 2. Excludes rows where:
#    - Species is Bombus barbutellus, B. rupestris, or B. vestalis
#    - Sequence is empty
#    - Gene assignment (from Primer) is NA
# 3. Assigns:
#    - "Gene" from Primer
#    - "Orientation" from Primer suffix ("F" = forward, "R" = reverse)
# 4. Adds columns:
#    - Gene
#    - Orientation
#    - fastaHeader (custom formatted FASTA header):
#      >ID|Sex|Orientation [organism=Species]
################################################################################

if [ -z "$1" ]; then
  echo "Usage: ./01_cleanDatabase.sh inputfile.csv"
  exit 1
fi

INPUT="$1"
OUTPUT="${INPUT%.csv}_filtered.csv"

tr -d '\r' < "$INPUT" | awk -F',' -v OFS=',' '
NR == 1 {
  print $1, $2, $3, $4, $5, $6, $7, $8, "Gene", "Orientation", "fastaHeader"
  next
}
{
  species = $2
  primer = $4
  sequence = $8

  if (species == "Bombus barbutellus" || species == "Bombus rupestris" || species == "Bombus vestalis")
    next

  if (sequence == "") next

  # Assign gene
  gene = "NA"
  if (primer ~ /^Hymeno[[:space:]]*[12]?$/ || primer ~ /^Lep[[:space:]]*[FR]?$/)
    gene = "COX1"
  else if (primer ~ /^16[[:space:]]*S/)
    gene = "16S"

  if (gene == "NA") next

  # Assign orientation
  orientation = "NA"
  if (primer ~ /F$/ || primer ~ /WB$/ || primer ~ /1$/)
    orientation = "F"
  else if (primer ~ /R$/ || primer ~ /IR$/ || primer ~ /2$/)
    orientation = "R"

  # fasta header:
  # >Gene|ID|Species|Sex|Primer|Orientation|Sequencing|Code|Nanodrop [organism=Species]
  fasta = ">" $1 "|" $7 "|" orientation " [organism="$2"]"

  print $1, $2, $3, $4, $5, $6, $7, $8, gene, orientation, fasta
}
' > "$OUTPUT"

echo "Filtered file saved to: $OUTPUT"
