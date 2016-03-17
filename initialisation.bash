#!/usr/bin/env bash

set -u

# Check for the right number of arguments. Assign them to variables.
NumArgsExpected=3
if [ "$#" -ne "$NumArgsExpected" ]; then
  echo "$#" 'arguments specified;' "$NumArgsExpected" 'expected. Quitting' >&2
  exit 1
fi
ConfigFile="$1"
RefAlignment="$2"
OutDir="$3"

# Source required code & check files exist
ThisDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source "$ThisDir"/'shiver_funcs.bash'
CheckFilesExist "$ConfigFile" "$RefAlignment"
source "$ConfigFile"

# If OutDir does not exist, try to create it.
if [ ! -d "$OutDir" ]; then
  mkdir "$OutDir"
fi || { echo 'Unable to create the specified output directory. Quitting.' >&2 ;\
exit 1; }

# Remove a trailing slash, if present.
OutDir=$(cd "$OutDir"; pwd)

# Check OutDir is empty
if ! find "$OutDir"/ -maxdepth 0 -empty | read v; then
  echo "$OutDir" 'is not empty. Delete or move its contents. Quitting.' >&2
  exit 1;
fi

# Some files we'll create
NewRefAlignment="$OutDir"/'ExistingRefAlignment.fasta'
RefList="$OutDir"/'ExistingRefNamesSorted.txt'
UngappedRefs="$OutDir"/'ExistingRefsUngapped.fasta'
database="$OutDir"/'ExistingRefsBlastDatabase'

# TODO: strip pure-gap columns - global alignment reconstruction depends on it.
cp "$RefAlignment" "$NewRefAlignment"

# List all names in the reference alignment
awk '/^>/ {print substr($1,2)}' "$NewRefAlignment" | sort > "$RefList"

# Check that RefAlignment has some sequences, and that their IDs are unique.
NumRefs=$(wc -l "$RefList" | awk '{print $1}')
if [[ $NumRefs -eq 0 ]]; then
  echo "$RefAlignment contains no sequences. Quitting." >&2
  exit 1;
fi
NumUniqueIDs=$(uniq "$RefList" | wc -l)
if [[ $NumUniqueIDs -ne $NumRefs ]]; then
  echo "$RefAlignment contains some identically named sequences. Rename"\
  'these and try again. Quitting.' >&2
  exit 1;
fi

# Ungap RefAlignment
"$Code_UngapFasta" "$NewRefAlignment" > "$UngappedRefs" || \
{ echo 'Problem ungapping' "$RefAlignment"'. Quitting.' >&2 ; exit 1; }

# Check that OutDir does not have whitespace in it
if [[ "$OutDir" =~ ( |\') ]]; then
  echo 'Unfortunately, blast cannot handle whitespace in paths.'\
  '(Stupid, I know.) This script therefore cannot execute the command' >&2
  echo "$BlastDBcommand" -dbtype nucl -in "$UngappedRefs" -input_type fasta \
  -out "$database" >&2
  echo "You'll have to use a directory without whitespace, or else move" \
  "$UngappedRefs to another directory, run the above command, and move its" \
  "output back to $OutDir. Sorry about that. Quitting." >&2
  exit 1;
fi

# Create the blast database
"$BlastDBcommand" -dbtype nucl -in "$UngappedRefs" -input_type fasta -out \
"$database" || \
{ echo 'Problem creating a blast database out of' \
"$OutDir"/'ExistingRefsUngapped.fasta. Quitting.' >&2 ; exit 1; }
