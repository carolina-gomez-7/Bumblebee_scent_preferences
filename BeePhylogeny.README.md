____
**DATE:**  20250613 
**AUTHOR:** CGR & MZF
_____
**GOAL :**
To build a phylogeny of Bumblebee species using sanger sequencing data from 2 different genes. 
____

# PATH
```
/Users/mzuluagaf/Documents/PhD_obsidian/04_Personal/Carolina/Bee_phylogeny
```

# 1) Filter, annotate, and prepare bumblebee Sanger sequence data
Since the data from all bumblebees is in a single database `Bumblebee_Durham_all.csv` , we first have to curate a few elements in it: 
1. species that will be included
2. Identify wich genes correspond to wich sequences.
3. identify the orientation of the sequences. 
4. Generate a header for the future FASTA files that will be generated. 
These procedures are done with the script `01_cleanDatabase.sh` that is stored in the scripts folder. 
The way to call it is by using the following comand:
```bash
bash ../scripts/01_cleanDatabase.sh Bumblebee_Durham_all.csv
```
It will generate an output file with the same name as the input but followed by the `_filtered.csv` sufix.

## 01_cleanDatabase.sh
```bash title:01_cleanDatabase.sh
#!/bin/bash

################################################################################
# Script:     01_cleanDatabase.sh
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
  fasta = ">" $1 "|" $7 "|" orientation [organism="$2"]"

  print $1, $2, $3, $4, $5, $6, $7, $8, gene, orientation, fasta
}
' > "$OUTPUT"

echo "Filtered file saved to: $OUTPUT"

```
# 2) Generate FASTA file(s) from filtered bumblebee Sanger CSV database
Using the fastaHeader column and the sequences, we can generate a FASTA file.
The following script uses the `fastaHeader` column and the `sequence` column for generating these .fasta files.
By default, the script generates a single fasta file containing all the sequences available in the database.
Custom FASTA files with subsets of the sequences can also vbe retrieved. Thes subsets can be filtered using the parameters (columns) from the database. 
In order to apply a filter, it is necessary to specify it by using a `flag`.
A list of flags is avilable at the beggining of the script and can be called by simply calling the script `bash 02_generateFasta.sh`.
A suggested implementation of the script where the database is filtered by gene and position is detailed below. Foru files should be generated form this implementation.
```bash 
# generate fasta from the database filtering by gene and orientation
bash ../scripts/02_generateFasta.sh Bumblebee_Durham_all_filtered.csv -g -or
```
## 02_generateFasta.sh
``` bash title:02_generateFasta.sh
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

```

# 03) Fix the orientation and get one single file per gene
For each gene, there are two primers (forward and reverse) hence sequencing have read the PCR amplicons in two two opposite orientations.
Some aligner algorithms do not reverse the sequences by default, therefore, it is necessary to convert one of the orientations to its reverse complement and only then, to concatenate it with the forward sequences of that given gene.    
This step is done manually using the online tool [seqret](https://www.ebi.ac.uk/jdispatcher/sfc/emboss_seqret) for the reverse sequences. It should also work for the Forward sequences. The orientation used here was selected for convenience only.
```sh
# these are the parameters used in the webpage (for reference)
# 1) Load the reverse strands of both genes and save them (manually) with a sufix 
singularity exec singularity/emboss:6.6.0 /emboss/bin/seqret -auto -stdout -sequence emboss_seqret-I20250614-074045-0122-69755028-p1m.sequence -snucleotide1 -sformat1 pearson -osformat2 fasta -feature -ofname2 emboss_seqret-I20250614-074045-0122-69755028-p1m.gff -sreverse1

singularity exec singularity/emboss:6.6.0 /emboss/bin/seqret -auto -stdout -sequence emboss_seqret-I20250614-075540-0639-3796619-p1m.sequence -snucleotide1 -sformat1 pearson -osformat2 fasta -feature -ofname2 emboss_seqret-I20250614-075540-0639-3796619-p1m.gff -sreverse1
```

Once the reverse sequences are reversed and only then, these reverSED sequences will be concatenated with the forward sequences of their corresponding gene. 
Save the reversed sequence files in the same working folder using the suffix `_reversed.fasta`
Use the following hardcoded script to concatenate the sequences. 
```bash
bash 03_fixOrientation_and_linearizeFasta.sh
```

## 03_mergeGenes_and_linearizeFasta.sh
```bash title:03_fixOrientation_and_linearizeFasta.sh
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
```

# 04) Align each gene (using Clustal Omega - online)
For each gene, alignement of all the sequences will be performed using the ClustalOmega algorithm. This algorithm is efficient for this amount of sequences.
Follow these steps:
1) Use the **EMBL-EBI web interface**: [Clustal Omega online](https://www.ebi.ac.uk/jdispatcher/msa/clustalo?outfmt=fa)
2) Use pearson/FASTA as the output format
3) Download using the option `Alignment in FASTA format`, since this keeps the names in a correct shape

# 05) Trim sequences to obtain common coordinates
Now that the sequences for each gene are aligned, trimming of the sequence edges can be performed in order to have common coordinates 
Comon coordinates that only retain phylogenetically informative regions will allow for colapsing each species into a consensus sequence in the folllowing steps. 
Start by downloading the software Gblocks (manually) and save it in the scripts folder. 
Then use it to trim the sequences of each gene.
```bash
# 1) download (manually) Gblocks from https://kbase.us/applist/apps/kb_gblocks/run_Gblocks/release
# 1.1) Move the file to the working folder (for each gene)

# 2) Trim sequences using Gblocks
/Users/mzuluagaf/Documents/PhD_obsidian/04_Personal/Carolina/Bee_phylogeny/scripts/Gblocks_0.91b/Gblocks COX1_clustal.aln-fasta -t=d -b5=h -b4=5
/Users/mzuluagaf/Documents/PhD_obsidian/04_Personal/Carolina/Bee_phylogeny/scripts/Gblocks_0.91b/Gblocks 16S_clustal.aln-fasta -t=d -b5=h -b4=5

```

# 6) Clean the headers from all sequences 
Now that all sequences share the same coordinates (per gene), we should clean the fasta headers to only retain the species information. This will simplify the process of sequence collapsing into consensus sequences.  
For each output file form the previous script, clean the headers using the follosing script:
```bash
bash ./04_cleanHeaders.sh input_file.gb
```

## 06_cleanHeaders.sh
```bash title:06_cleanHeaders.sh
#!/bin/bash
##################################################################################
# Script:     06_cleanHeaders.sh
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

awk '/^>/ {
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

```

# 5) Concatenate Genes
Since each sequence in each gene has the same coordinates, we can use AMAS.py sumarize them into single consensus sequences. Use the following lines directly into the terminal. 
No bash script is provided this time since some debugging may be neccesary for installing the software. 
```bash
# 1) Install AMAS <- Necessary for concatenting both genes
# Download and place the working folder
brew install git
brew install python
git clone https://github.com/marekborowiec/AMAS.git

# 5) concatenate genes
python3 /Users/mzuluagaf/Documents/PhD_obsidian/04_Personal/Carolina/Bee_phylogeny/data/AMAS/amas/AMAS.py concat -f fasta -d dna -i COX1_clustal_clean.aln-fasta-gb 16S_clustal_clean.aln-fasta-gb -p partitions.txt -t concatenated.fas

sed -i '' 's/^/DNA, /' partitions.txt

```

> [!summary] Process so far:
>  # What You Now Have After Concatenation
> Each row in your concatenated alignment represents a species, and the sequence for that row is made of:
> - The aligned, trimmed sequence of gene A,
> - directly followed by
> - The aligned, trimmed sequence of gene B.
> This creates a composite (multi-locus) representation for each species.
> 
> # Conceptually:
> Each species is represented by a pseudo-sequence that is the concatenation of orthologous gene regions from that species.
> Even though the sequences may have come from different individuals within a species, this is standard practice in phylogenetics when building species-level trees using Sanger data, especially when:
> - Individuals are confidently assigned to species (which you’ve done using BLAST + morphology).
> - Genes are orthologous and alignable across taxa.
> - There is no intention to infer within-species variation or population structure.
> 
> # Summary: What You're Modeling
> You're effectively building a species tree based on concatenated, orthologous gene sequences.
> This assumes that within-species variation is negligible or unimportant for your goal.
> Each species’ composite sequence captures evolutionary signal from multiple loci, improving phylogenetic resolution.



# 07_buildPhylogeny
```bash
# Istall some programs:
brew install eigen
# check the location of eigen
brew info eigen
# Intall boost
brew install boost
# 1) Install IQ-tree using git
git clone --recurse-submodules https://github.com/iqtree/iqtree2.git
cd iqtree2
mkdir build && cd build
cmake \
  -DEIGEN3_INCLUDE_DIR=/opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3 \
  -DBOOST_ROOT=/opt/homebrew/opt/boost \
  ..
make -j4

# download manually from https://iqtree.github.io then place the folder in wuyour data directory and use its paht to run the program 
/Users/mzuluagaf/Documents/PhD_obsidian/04_Personal/Carolina/Bee_phylogeny/data/iqtree-3.0.1-macOS/bin/iqtree3 -h

# Run IQ-Tree
/Users/mzuluagaf/Documents/PhD_obsidian/04_Personal/Carolina/Bee_phylogeny/data/iqtree-3.0.1-macOS/bin/iqtree3 -s concatenated.fas -p partitions.txt -m MFP+MERGE -B 1000 -T AUTO
```

```bash 
# Check for best model
> cat partitions.txt.best_scheme.nex
#nexus
begin sets;
  charset p1_COX1_clustal_clean = 1-544;
  charset p2_16S_clustal_clean = 545-948;
  charpartition mymodels =
    TIM2+F+I: p1_COX1_clustal_clean,
    TIM3+F+I: p2_16S_clustal_clean;
end;

```

> [!NOTE]- INTERPRETATION
> 
> partitions.txt.best_scheme.nex — Explained
> 
> #nexus
> begin sets;
>   charset p1_COX1_clustal_clean = 1-544;
>   charset p2_16S_clustal_clean = 545-948;
> ✅ This defines two partitions:
> 
> p1_COX1_clustal_clean: sites 1–544 in the concatenated alignment
> p2_16S_clustal_clean: sites 545–948
> So your concatenated alignment has two gene regions, exactly as expected.
> 
>   charpartition mymodels =
>     TIM2+F+I: p1_COX1_clustal_clean,
>     TIM3+F+I: p2_16S_clustal_clean;
> end;
> ✅ This defines a partitioning scheme named mymodels with:
> 
> TIM2+F+I applied to COX1
> TIM3+F+I applied to 16S
> 🧠 Model meanings:
> TIM2 / TIM3: Transitional models (more flexible than HKY, simpler than GTR)
> +F: Uses observed empirical base frequencies
> +I: Includes a proportion of invariable sites
> 🟨 These models were automatically selected by IQ-TREE as best-fitting per partition using BIC (Bayesian Information Criterion).
> 
> ✅ Summary of Interpretation
> 
> IQ-TREE determined that:
> COX1 evolves best under TIM2+F+I from sites 1–544
> 16S evolves best under TIM3+F+I from sites 545–948
> and did not merge the partitions, treating them separately.
> 🧾 What to Report in Your Methods
> 
> You can write something like:
> 
> "We used IQ-TREE v3.0.1 to perform maximum likelihood phylogenetic inference under a partitioned model. Substitution models were selected independently for each gene region using ModelFinder and the Bayesian Information Criterion (BIC). The best-fit models were TIM2+F+I for the COX1 partition (sites 1–544) and TIM3+F+I for the 16S partition (sites 545–948). Bootstrap support was assessed using 1000 ultrafast bootstrap replicates."
>
>**You do _not_ need to re-run IQ-TREE again** if you used:
>```
>iqtree3 -s concatenated.fas -p partitions.txt -m MFP+MERGE -B 1000 -T AUTO
>```
>**Because:**
> `-m MFP+MERGE` already:
> - Tested many models for each partition,
> - Selected the best-fit model per partition using BIC,
> - Optionally merged partitions if warranted,
> - Then built the tree using those best-fit models.

# RESULTS:
## Rooted
![[Pasted image 20250615002939.png]]

## Unrooted

![[Pasted image 20250615004424.png]]