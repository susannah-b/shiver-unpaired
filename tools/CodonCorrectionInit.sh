#!/usr/bin/env bash

set -u
set -o pipefail

source "$HOME/shiver/config.sh" # adapt for flexible shiver location

# This script is to be run once prior to the use of CodonCorrection.sh.
# Could eventually be integrated as part of the shiver init.
# Give specific description for gene coords
UsageInstructions="Arguments for this script:
(1) A fasta file of annotatable references
(2) The output folder for the BLAST database and init steps
(3) A file of gene coordinates
"

# Check correct arguments
NumArgsExpected=3
if [ "$#" -ne "$NumArgsExpected" ]; then
  echo "$UsageInstructions"
  echo "$#" 'arguments specified;' "$NumArgsExpected" 'expected. Quitting.' >&2
  exit 1
fi

ReferenceFasta="$1"
OutputDirInit="$2"
GeneCoordInfo="$3"

genes=("GAG" "POL" "VIF" "VPR" "VPU" "ENV" "NEF") # This is currently defined in both the init and CC.sh

# Variables # Set up for $shiver
PythonFuncs="$HOME/shiver/tools/CC_Python_funcs.py" # for full implementation add this to shiver init so i can have one variable for CC.sh and CC Init

# Check necessary files exist TODO

# Check the output directory and cd
if [[ -d "$OutputDirInit" ]]; then
  echo "$OutputDirInit exists already; quitting to prevent overwriting." >&2
  exit 1
fi
mkdir -p "$OutputDirInit" && cd "$OutputDirInit" ||
{ echo "Could not mkdir then cd to $OutputDirInit. Quitting." >&2 ; exit 1 ; }


### Check the reference file is a fasta
if [[ -f "$ReferenceFasta" ]]; then
# Check reference file is a .fasta
  if [[ "$ReferenceFasta" != *.fasta ]]; then 
    echo "Reference file $ReferenceFasta is not a fasta file. Quitting." >&2
    exit 1
  fi

  # Check the reference file contains sequences
  Ref_SeqNumber=$(grep '^>' "$ReferenceFasta" | wc -l | awk '{$1=$1};1')
  if [[ "$Ref_SeqNumber" == 0 ]]; then
    echo "Reference file $ReferenceFasta contains no sequences. Quitting." >&2
    exit 1
  fi
else
  echo "Reference file $ReferenceFasta does not exist. Check specified filepath. Quitting." >&2
  exit 1
fi

# Count sequence number in reference file
Ref_SeqNumber=$(grep '^>' "$ReferenceFasta" | wc -l | awk '{$1=$1};1')
if [[ "$Ref_SeqNumber" == 0 ]]; then
  echo "Reference file $ReferenceFasta contains no sequences. Quitting." >&2
  exit 1
fi

# Check GeneCoordInfo
if [[ "$GeneCoordInfo" != *.txt ]]; then 
  echo "Gene coordinates file in unexpected format. Please provide as a .txt file. Quitting." >&2
  exit 1
fi

# Check GeneCoordInfo header is as expected
# Not strictly necessary for function but ensures the gene info is formatted in the correct way. But maybe overly strict.
GC_FirstLine=$(head -n 1 "$GeneCoordInfo")
GeneCoord_Header="Sequence_name, gag_start, gag_end, pol_start, pol_end, vif_start, vif_end, vpr_start, \
vpr_end, vpu_start, vpu_end, env_start, env_end, nef_start, nef_end"
if [[ "$GC_FirstLine" != *"$GeneCoord_Header"* ]]; then
  echo -e "Expected the gene coordinates file to have the header: '$GeneCoord_Header'\nPlease supply the gene coordinate \
  data in the correct format. Quitting" >&2
  exit 1
fi

# Extract genes for each gene
$python2 "$PythonFuncs" MakeReferenceDatabase --InitDir "$OutputDirInit" --GeneCoordInfo "$GeneCoordInfo" --GenomeFile "$ReferenceFasta" || { echo "CodonCorrectionInit.sh was unable to extract the reference gene sequences. Quitting." >&2; exit 1; }

# Copy reference file and Gene Coords to working directory
Ref_Copied="$OutputDirInit"/CC_References.fasta
cp "$ReferenceFasta" "$Ref_Copied" || { echo "Failed to copy $ReferenceFasta to the output directory. Quitting." >&2; exit 1; }
GeneCoords_Copied="$OutputDirInit"/CC_Coords.fasta
cp "$GeneCoordInfo" "$GeneCoords_Copied" || { echo "Failed to copy $GeneCoordInfo to the output directory. Quitting." >&2; exit 1; }

# Create a whole genome database to BLASTn to
for gene in "${genes[@]}"; do
  BLASTnDB="$OutputDirInit"/BLASTnDB_${gene}
  if [ ! -d "$BLASTnDB" ]; then
    mkdir "$BLASTnDB"
  fi || { echo "Unable to create database folder. Quitting."  >&2 ; exit 1 ; }

  # Make the BLASTn database
  GeneFile="$OutputDirInit/ReferenceGenes_${gene}.fasta"
  makeblastdb -dbtype nucl -in "$GeneFile" -input_type fasta -out "$BLASTnDB"/"${gene}" || { echo "Unable to create \
  a blast database in $BLASTnDB. Quitting." ; exit 1; }
done

# Extract HXB2 genes needed in CodonCorrection.sh to extract sample genes
B.FR.83.HXB2_LAI_IIIB_BRU.K03455, 790, 2292, 2253, 5096, 5041, 5619, 5559, 5850, 6062, 6310, 6225, 8795, 8797, 9168, 

# Can add to this init depending on what data is included with shiver. Currently I have a separate pipeline which takes Nick's all-cds-info file and extracts the accession numbers which are 
# then used to download the full sequences from LANL. Then it uses the supplied gene fragments with the full genomes to extract gene coordinates (with a separate step for extracting Pol from Gag-Pol)
# Current CodonCorrection.sh usage assumes the presence of a whole genome reference file and a list of gene coordinates in the exact correct format. 
