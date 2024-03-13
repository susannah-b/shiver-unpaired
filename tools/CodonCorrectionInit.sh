#!/usr/bin/env bash

set -u
set -o pipefail

# This script is to be run once prior to the use of CodonCorrection.sh.
# Could eventually be integrated as part of the shiver init.

# For testing purposes; copy command below (specific to my machine):
# /Users/s.blundell/shiver-unpaired/tools/CodonCorrectionInit.sh ~/shiver/CodonCorrectionReferences/HIV1_ALL_2020_genome_DNA_no_ambiguity_only_annotated_aligned_no_gap_columns_aligned.fasta

UsageInstructions="Arguments for this script:
(1) A fasta file of annotatable references
(2) The output folder for the BLAST database and init steps
"

NumArgsExpected=2
if [ "$#" -ne "$NumArgsExpected" ]; then
  echo "$UsageInstructions"
  echo "$#" 'arguments specified;' "$NumArgsExpected" 'expected. Quitting.' >&2
  exit 1
fi
ReferenceFasta="$1"
OutputDirInit="$2" # Not fully implemented this yet - should be read by CC.sh

# Check necessary files exist TODO
# cd to $OutputDir
if [[ -d "$OutputDirInit" ]]; then
  echo "$OutputDirInit exists already; quitting to prevent overwriting." >&2
  exit 1
fi
mkdir -p "$OutputDirInit" && cd "$OutputDirInit" ||
{ echo "Could not mkdir then cd to $OutputDirInit. Quitting." >&2 ; exit 1 ; }

if [[ -f "$ReferenceFasta" ]]; then
# Check reference file is a .fasta
  if [[ "$ReferenceFasta" != *.fasta ]]; then 
    echo "Reference file $ReferenceFasta is not a fasta file. Quitting." >&2
    exit 1
  fi

  ### Check the reference file contains sequences
  # Count sequence number in reference and sample
  Ref_SeqNumber=$(grep '^>' "$ReferenceFasta" | wc -l | awk '{$1=$1};1')
  if [[ "$Ref_SeqNumber" == 0 ]]; then
    echo "Reference file $ReferenceFasta contains no sequences. Quitting." >&2
    exit 1
  fi
else
  echo "Reference file $ReferenceFasta does not exist. Check specified filepath. Quitting." >&2
  exit 1
fi

BLASTnDB=""$OutputDirInit"/BLASTnDB"
if [ ! -d "$BLASTnDB" ]; then
  mkdir "$BLASTnDB"
fi || { echo "Unable to create database folder. Quitting."  >&2 ; exit 1 ; }

# Copy reference file to working directory
# figure out non-specific shiver path
Ref_Copied=""$OutputDirInit"/CC_References.fasta"
cp "$ReferenceFasta" "$Ref_Copied" || { echo "Failed to copy $ReferenceFasta to the output"\
  "directory. Quitting." >&2; exit 1; }

### Create a whole genome database to BLASTn to
# should really call the reference file as an arg but set that up later
# filepath to db is not generic yet
makeblastdb -dbtype nucl -in "$ReferenceFasta" -input_type fasta -out "$BLASTnDB"/CC_whole_genome || { echo "Unable to create \
a blast database in "$BLASTnDB". Quitting." ; exit 1; }