#!/usr/bin/env bash

set -u
set -o pipefail

# Find shiver files
ToolsDir="$(cd "$(dirname "${BASH_SOURCE[0]}" )" && pwd)"
shiver_bin="$(dirname "$ToolsDir")" # /path/to/shiver/bin
source "$shiver_bin"/'config.sh'
source "$shiver_bin"/'shiver_funcs.sh'
PythonFuncs="$shiver_bin"/'tools/CC_Python_funcs.py'

# This script is to be run once prior to the use of CodonCorrection.sh. The current version of the init assumes the presence of a whole genome 
# reference file and a list of gene coordinates in the exact correct format - header included. 

# Give specific description for gene coords
GeneCoord_Header="Sequence_name, gag_start, gag_end, pol_start, pol_end, vif_start, vif_end, vpr_start, \
vpr_end, vpu_start, vpu_end, env_start, env_end, nef_start, nef_end"

UsageInstructions="Arguments for this script:
(1) A fasta file of aligned annotatable references
(2) The output folder for the BLAST database and init steps
(3) A file of gene coordinates as a txt file, with the header: '$GeneCoord_Header'
The gene coordinates for each sequence should be listed per line in the same order as the above header, separated by commas.
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

CheckFilesExist "$ReferenceFasta" "$GeneCoordInfo"

# Check the output directory and cd
if [[ -d "$OutputDirInit" ]]; then
  echo "$OutputDirInit exists already; quitting to prevent overwriting." >&2
  exit 1
fi
mkdir -p "$OutputDirInit" && cd "$OutputDirInit" ||
{ echo "Could not mkdir then cd to $OutputDirInit. Quitting." >&2 ; exit 1 ; }

# Check the reference file contains sequences
if ! grep -e '^>' "$ReferenceFasta" > /dev/null; then
  echo "Reference file $ReferenceFasta contains no sequences. Quitting." >&2
  exit 1
fi

# Check GeneCoordInfo header is as expected
  # Not strictly necessary for function but ensures the gene info is formatted in the correct way.
GC_FirstLine=$(head -n 1 "$GeneCoordInfo")
if [[ "$GC_FirstLine" != *"$GeneCoord_Header"* ]]; then
  echo -e "Expected the gene coordinates file to have the header: '$GeneCoord_Header'\nPlease supply the gene coordinate \
  data in the correct format and order. If coordinates are already supplied in the correct order, update the header to match. Quitting" >&2
  exit 1
fi

# Extract each gene to make a BLAST database for each
$python "$PythonFuncs" MakeReferenceDatabase --InitDir "$OutputDirInit" --GeneCoordInfo "$GeneCoordInfo" --GenomeFile "$ReferenceFasta" \
|| { echo "Failed to extract the reference gene sequences. Quitting." >&2; exit 1; }

# Copy reference file and gene coordinates file to working directory
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
  "$BlastDBcommand" -dbtype nucl -in "$GeneFile" -input_type fasta -out "$BLASTnDB"/"${gene}" -logfile "$BLASTnDB/$gene.log" || { echo "Unable to create \
  a blast database in $BLASTnDB. Check that the shiver config has the correct file path for 'BLASTDBcommand'. Quitting." ; exit 1; }
done
