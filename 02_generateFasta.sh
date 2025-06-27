#!/bin/bash

################################################################################
# Date:      20250616    
# Script:     02_generateFasta.sh
# Author:     Carolina Gomez Ramirez & Maximiliano Zuluaga Forero
# Purpose:    Generate FASTA file(s) from filtered bumblebee Sanger CSV database
# Usage:      ./02_generateFasta.sh input_filtered.csv [flags]
#
# FLAGS:
#   -g     Gene
#   -id    ID
#   -sp    Species
#   -s     Sex
#   -p     Primer
#   -or    Orientation (F or R)
#   -seq   Sequencing code
#   -nd    Nanodrop filter (e.g., =62.5 or >50 or <80) <- NOTE: This function
#                                                         hasn't been test fully
#
# Behavior:
#   - Without flags: full FASTA of all entries
#   - With flags: expands to one FASTA per group combination
#   - Multiple flags accepted, will filter in order of the flags
###############################################################################%

if [ -z "$1" ]; then
  echo -e "\n❗ No input file provided ❗ \n"
  echo -e "Usage: ./02_generateFasta.sh input_filtered.csv [flags]"
  echo -e "bash   ^^^^ script.sh ^^^^^   ^^^ input.csv ^^^ [flags]   <<---------- EXAMPLE \n"
  awk 'NR==1,/^###############################################################################%$/' "$0"
  exit 1
fi


INPUT="$1"
shift

# Initialize filter variables
GENE="" ID="" SPECIES="" SEX="" PRIMER="" ORIENT="" SEQCODE="" NDOP=""
EXPAND=()

# Parse flags
while [[ $# -gt 0 ]]; do
  case "$1" in
    -g)    if [[ -n "$2" && "$2" != -* ]]; then GENE="$2"; shift; else EXPAND+=("Gene"); fi ;;
    -id)   if [[ -n "$2" && "$2" != -* ]]; then ID="$2"; shift; else EXPAND+=("ID"); fi ;;
    -sp)   if [[ -n "$2" && "$2" != -* ]]; then SPECIES="$2"; shift; else EXPAND+=("Species"); fi ;;
    -s)    if [[ -n "$2" && "$2" != -* ]]; then SEX="$2"; shift; else EXPAND+=("Sex"); fi ;;
    -p)    if [[ -n "$2" && "$2" != -* ]]; then PRIMER="$2"; shift; else EXPAND+=("Primer"); fi ;;
    -or)   if [[ -n "$2" && "$2" != -* ]]; then ORIENT="$2"; shift; else EXPAND+=("Orientation"); fi ;;
    -seq)  if [[ -n "$2" && "$2" != -* ]]; then SEQCODE="$2"; shift; else EXPAND+=("Sequencing code"); fi ;;
    -nd)   NDOP="$2"; shift ;;
    *)     echo "Unknown flag: $1"; exit 1 ;;
  esac
  shift
done

# No filters or expansions → dump full file
if [[ -z "$GENE$ID$SPECIES$SEX$PRIMER$ORIENT$SEQCODE$NDOP" && ${#EXPAND[@]} -eq 0 ]]; then
  awk -F',' '
    NR==1 { for (i=1;i<=NF;i++) h[$i]=i; next }
    { print $h["fastaHeader"] "\n" $h["Sequences"] }
  ' "$INPUT" > "${INPUT%.csv}.fasta"
  echo "✔ Wrote: ${INPUT%.csv}.fasta"
  exit 0
fi

# If any fixed filter is used
if [[ -n "$GENE$ID$SPECIES$SEX$PRIMER$ORIENT$SEQCODE$NDOP" ]]; then
  awk -F',' -v gene="$GENE" -v id="$ID" -v sp="$SPECIES" -v sex="$SEX" \
      -v primer="$PRIMER" -v orient="$ORIENT" -v seqcode="$SEQCODE" -v ndop="$NDOP" '
    BEGIN { IGNORECASE=1 }
    NR==1 { for (i=1;i<=NF;i++) h[$i]=i; next }
    {
      if (gene && $h["Gene"] != gene) next
      if (id && $h["ID"] != id) next
      if (sp && $h["Species"] != sp) next
      if (sex && $h["Sex"] != sex) next
      if (primer && $h["Primer"] != primer) next
      if (orient && $h["Orientation"] != orient) next
      if (seqcode && $h["Sequencing code"] != seqcode) next

      if (ndop != "") {
        n = $h["Nanodrop concentration"] + 0
        op = substr(ndop, 1, 1)
        val = substr(ndop, 2)
        if (op == ">" && !(n > val)) next
        if (op == "<" && !(n < val)) next
        if (op == "=" && !(n == val)) next
      }

      print $h["fastaHeader"] "\n" $h["Sequences"]
    }
  ' "$INPUT" > "filtered_${INPUT%.csv}.fasta"

  echo "✔ Wrote: filtered_${INPUT%.csv}.fasta"
  exit 0
fi

# Multi-group expansion mode
awk -F',' -v OFS=',' -v cols="${EXPAND[*]}" -v inputname="$(basename "$INPUT" .csv)" '
BEGIN {
  ncols = split(cols, colArr, " ")
}
NR==1 {
  for (i=1; i<=NF; i++) {
    col = $i
    gsub(/^ +| +$/, "", col)  # Trim spaces
    h[col] = i
  }
  next
}
{
  key = ""
  for (i = 1; i <= ncols; i++) {
    colname = colArr[i]
    colidx = h[colname]
    if (colidx == 0) next  # skip if not found
    val = $colidx
    gsub(/[ \/]/, "_", val)
    key = (i == 1 ? val : key "_" val)
  }

  file = key "_" inputname ".fasta"
  lines[file] = lines[file] $h["fastaHeader"] "\n" $h["Sequences"] "\n"
}
END {
  for (f in lines) {
    print lines[f] > f
    close(f)
    print "✔ Wrote: " f
  }
}
' "$INPUT"


