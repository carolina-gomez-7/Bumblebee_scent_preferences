#!/bin/bash

################################################################################
#
# Date:      20250616    
# Script:     03_fixOrientation_and_linearizeFasta.sh
# Author:     Carolina Gomez Ramirez & Maximiliano Zuluaga Forero
# Purpose:    make files per gene with a single orientation (forward)
# Usage:      ./03_mergeGenes_and_linearizeFasta.sh
#
################################################################################

# 1) Go to seqret webpage >> https://www.ebi.ac.uk/jdispatcher/sfc/emboss_seqret << and convert the reverse 
# sequences into their reverse orientation, then download them and save the files in the same working folder using the suffix `_reversed.fasta`

# 2) Merge the sequences with the correct orientation from the same gene in a single file
cat 16S_F_Bumblebee_Durham_all_filtered.fasta 16S_R_Bumblebee_Durham_all_filtered_reversed.fasta > 16S_bumblebee_durham_all_filtered.fasta
cat COX1_F_Bumblebee_Durham_all_filtered.fasta COX1_R_Bumblebee_Durham_all_filtered_reversed.fasta > COX1_bumblebee_durham_all_filtered.fasta 

# 3) Make fastas in single line files
seqtk seq 16S_bumblebee_durham_all_filtered.fasta > 16S_bumblebee_durham_all_filtered_oriented.fasta
seqtk seq COX1_bumblebee_durham_all_filtered.fasta > COX1_bumblebee_durham_all_filtered_oriented.fasta
