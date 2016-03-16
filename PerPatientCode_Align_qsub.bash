#!/usr/bin/env bash

#PBS -l walltime=0:04:59
#PBS -l select=1:ncpus=1:mem=6GB
#PBS -q pqeelab
#PBS -J 1-2
PBS_ARRAY_INDEX=2

set -u

# TODO: make sure ~/PerPatientCode_Align.bash is up to date each time

# A file containing the paths to all of the fasta files raw out of IVA: 
ListOfInputFiles="$HOME/JobInputs/FileNames_SA_IVAraw_SimpleSeqs.txt"
# for BEEHIVE: RefAlignment="$HOME/JobInputs/HIV1_COM_2012_genome_DNA.fasta"
RefAlignment="$HOME/JobInputs/HIV1_COM_2012_genome_DNA_NoLTR.fasta"
OutputDir="$WORK/test"

module load python/2.7.3
module load blast+/2.2.30
module load mafft/7
module load R/3.1.2

# Check that certain input files exist
function CheckFilesExist {
  for argument in "$@"; do
    if [ ! -f "$argument" ]; then
      echo "$argument" 'does not exist. Quitting.' >&2
      exit 1
    fi
  done
}
CheckFilesExist "$ListOfInputFiles" "$RefAlignment"

# Remove trailing slashes from the user-specified directories, if present.
if [ ! -d "$OutputDir" ]; then
  mkdir "$OutputDir"
fi
OutputDir=$(cd "$OutputDir"; pwd)

# Each line in the list of input files is an input file, so count the number of
# lines.
NumInputFiles=$(wc -l "$ListOfInputFiles" | awk '{print $1}')

# Check that the current array index is not larger than the number of input
# files.
if [ "$PBS_ARRAY_INDEX" -gt "$NumInputFiles" ]; then
  echo 'Error: job array index' "$PBS_ARRAY_INDEX" 'is greater than the number'\
  'of lines in' "$ListOfInputFiles"'. Quitting.' >&2
  exit 1
fi

# Find the input file corresponding to this array index. Check it exists. Copy
# it to the current directory.
ContigFile=$(sed -n "$PBS_ARRAY_INDEX"'p' "$ListOfInputFiles")
CheckFilesExist "$ContigFile"

# Set the contig file basename to be the file name minus extension and path.
# We'll use this for naming output files, so delete any similarly named files in
# the current directory. Then copy the contig file to the current directory.
ContigBaseName=$(basename "${ContigFile%.fasta}")
rm "$ContigBaseName"*
cp "$ContigFile" .
ContigFile=$(basename "$ContigFile")

# Info:
echo 'This is array job' "$PBS_ARRAY_INDEX" 'working with input file '\
"$ContigFile"

# The main program:
~/PerPatientCode_Align.bash "$ContigFile" "$RefAlignment" "$ContigBaseName" || \
{ echo 'Problem running PerPatientCode_Align.bash. Quitting.' >&2 ; exit 1 ; }

# Move output to the output directory
mv "$ContigBaseName"* -t "$OutputDir"

