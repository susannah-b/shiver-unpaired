#!/usr/bin/env bash

set -u

# Some files we'll create
UngappedRefs="$OutDir"/'ExistingRefsUngapped.fasta'
database="$OutDir"/'ExistingRefsBlastDatabase'

# Some code we'll need
ThisDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ToolsDir="$ThisDir"/tools
Code_UngapFasta="$ToolsDir/UngapFasta.py"

# Check for the right number of arguments. Assign them to variables.
NumArgsExpected=3
if [ "$#" -ne "$NumArgsExpected" ]; then
  echo "$#" 'arguments specified;' "$NumArgsExpected" 'expected. Quitting' >&2
  exit 1
fi
ConfigFile="$1"
RefAlignment="$2"
OutDir="$3"

# Check the ConfigFile exists. Source it.
if [ ! -f "$ConfigFile" ]; then
  echo "$ConfigFile" 'does not exist. Quitting.' >&2
  exit 1;
fi
source "$ConfigFile"

# Check the RefAlignment exists
if [ ! -f "$RefAlignment" ]; then
  echo "$RefAlignment" 'does not exist. Quitting.' >&2
  exit 1;
fi

# Check that RefAlignment has some sequences, and that their IDs are unique.
NumRefs=$(awk '/^>/ {print substr($1,2)}' "$RefAlignment" | wc -l)
if [[ $NumRefs -eq 0 ]]; then
  echo "$RefAlignment" 'contains no sequences. Quitting.' >&2
  exit 1;
fi
NumUniqueIDs=$(awk '/^>/ {print substr($1,2)}' "$RefAlignment" | sort | uniq | \
wc -l)
if [[ $NumUniqueIDs -ne $NumRefs ]]; then
  echo "$RefAlignment" 'contains some identically named sequences. Rename'\
  'these and try again. Quitting.' >&2
  exit 1;
fi

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

cp "$RefAlignment" "$OutDir"/'ExistingRefAlignment.fasta'

# Ungap RefAlignment
"$Code_UngapFasta" "$OutDir"/'ExistingRefAlignment.fasta' > "$UngappedRefs" || \
{ echo 'Problem ungapping' "$RefAlignment"'. Quitting.' >&2 ; exit 1; }

# Check that OutDir does not have whitespace in it
if [[ "$OutDir" =~ ( |\') ]]; then
  echo 'Unfortunately, blast cannot handle whitespace in paths.'\
  '(Stupid, I know.) This script therefore cannot execute the command' >&2
  echo "$BlastDBcommand" -dbtype nucl -in "$UngappedRefs" -input_type fasta \
  -out "$database" >&2
  echo "You'll have to use a directory without whitespace, or else move" \
  "$UngappedRefs" "to another directory, run the above command, and move it's" \
  "output back to $OutDir. Sorry about that. Quitting." >&2
  exit 1;
fi

# Create the blast database
"$BlastDBcommand" -dbtype nucl -in "$UngappedRefs" -input_type fasta -out \
"$database" || \
{ echo 'Problem creating a blast database out of' \
"$OutDir"/'ExistingRefsUngapped.fasta. Quitting.' >&2 ; exit 1; }
