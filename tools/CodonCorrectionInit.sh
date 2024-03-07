#!/usr/bin/env bash

set -u
set -o pipefail

# This script is to be run once prior to the use of CodonCorrection.sh.
# Could eventually be integrated as part of the shiver init.

# For testing purposes; copy command below (specific to my machine):
# /Users/s.blundell/shiver-unpaired/tools/CodonCorrectionInit.sh ~/shiver/CodonCorrectionReferences/HIV1_ALL_2020_genome_DNA_no_ambiguity_only_annotated_aligned_no_gap_columns_aligned.fasta

UsageInstructions="Arguments for this script:
(1) A fasta file of annotatable references
"

NumArgsExpected=1
if [ "$#" -ne "$NumArgsExpected" ]; then
  echo "$UsageInstructions"
  echo "$#" 'arguments specified;' "$NumArgsExpected" 'expected. Quitting.' >&2
  exit 1
fi
ReferenceFasta="$1"

# Specific filepath but I think functional as shiver is assumed to be in ~/shiver?
ReferenceFolder="$HOME/shiver/CodonCorrectionReferences/"
BLASTnDB="$HOME/shiver/CodonCorrectionReferences/BLASTnDB"

# Check necessary files exist TODO

### Create folders to store the database
# For now this file is created, if including a set of references with shiver you would presumably already create this folder
# Need to check other files required to be in this folder - it might be unecessary if there are now only 1-2 files
if [ ! -d "$ReferenceFolder" ]; then
	mkdir "$ReferenceFolder"
fi || { echo "Unable to create reference folder. Quitting."  >&2 ; exit 1 ; }

if [ ! -d "$BLASTnDB" ]; then
  mkdir "$BLASTnDB"
fi || { echo "Unable to create database folder. Quitting."  >&2 ; exit 1 ; }

### Create a whole genome database to BLASTn to
# should really call the reference file as an arg but set that up later
# filepath to db is not generic yet
makeblastdb -dbtype nucl -in "$ReferenceFasta" -input_type fasta -out "$BLASTnDB"/CC_whole_genome || { echo "Unable to create \
a blast database in "$BLASTnDB". Quitting." ; exit 1; }