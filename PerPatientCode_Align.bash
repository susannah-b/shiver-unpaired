#!/usr/bin/env bash

# Arguments for this script:
# 2) a fasta file of contigs (output from processing the short reads with an
# assembly program),
# 3) A sample ID ('SID') used for naming the output from this script (a sensible
# choice might be the contig file name minus its path and extension).
# If this script completes successfully it will produce a .blast file (detailing
# the contigs that have blast hits). If the .blast file is not empty it will
# produce a fasta file of those contigs with hits and another fasta file of
# these contigs aligned to the refereces; if cutting and/or reversing of contigs
# is necessary, two more fasta files are produced - the cut/reversed contigs on
# their own and also aligned to references (i.e. there will be two files of the
# contigs on their own and two files of the contigs aligned to references).

set -u

# Check for the right number of arguments. Assign them to variables.
NumArgsExpected=3
if [ "$#" -ne "$NumArgsExpected" ]; then
  echo "$#" 'arguments specified;' "$NumArgsExpected" 'expected. Quitting' >&2
  exit 1
fi
InitDir="$1"
ContigFile="$2"
SID="$3"

BlastDatabase="$InitDir/ExistingRefsBlastDatabase"
RefAlignment="$InitDir/ExistingRefAlignment.fasta"

# Source required code & check files exist
ThisDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source "$ThisDir"/'shiver_funcs.bash'
CheckFilesExist "$ContigFile" "$RefAlignment"

# Check that there are some contigs, and their IDs are unique.
ContigNames=$(awk '/^>/ {print substr($1,2)}' "$ContigFile" | sort)
NumContigs=$(echo "$ContigNames" | wc -w)
if [[ $NumContigs -eq 0 ]]; then
  echo "$ContigFile contains no sequences. Quitting." >&2
  exit 1;
fi
NumUniqueIDs=$(printf '%s\n' $ContigNames | uniq | wc -l)
if [[ $NumUniqueIDs -ne $NumContigs ]]; then
  echo "$ContigFile contains some identically named sequences. Rename"\
  'these and try again. Quitting.' >&2
  exit 1;
fi

# The names for output files we'll produce.
BlastFile="$SID"'.blast'
RawContigFile="$SID"'_hiv.fasta'
CutContigFile="$SID"'_hiv_cut.fasta'
RawContigAlignment="$SID"'_raw_wRefs.fasta'
CutContigAlignment="$SID"'_cut_wRefs.fasta'
TempRawContigAlignment='temp_'"$SID"'_raw_wRefs_swap.fasta'
TempCutContigAlignment='temp_'"$SID"'_cut_wRefs_swap.fasta'

# Blast the contigs
blastn -query "$ContigFile" -db "$BlastDatabase" -out "$BlastFile" \
-max_target_seqs 1 -outfmt \
'10 qacc sacc sseqid evalue pident qstart qend sstart send' || \
{ echo 'Problem blasting' "$ContigFile"'. Quitting.' >&2 ; exit 1 ; }

# If there are no blast hits, nothing needs doing. Exit.
NumBlastHits=$(wc -l "$BlastFile" | awk '{print $1}')
if [ "$NumBlastHits" -eq 0 ]; then
  echo "No contig in $ContigFile has a blast hit (i.e. this is presumably pure"\
  "contamination). Quitting."
  exit 0
fi

# Extract those contigs that have a blast hit...
HIVcontigNames=$(awk -F, '{print $1}' "$BlastFile" | sort | uniq)
"$Code_FindSeqsInFasta" "$ContigFile" $HIVcontigNames > "$RawContigFile" || \
{ echo 'Problem extracting the HIV contigs. Quitting.' >&2 ; exit 1 ; }

# ...and align them to the refs.
mafft --quiet --add "$RawContigFile" "$RefAlignment" \
> "$TempRawContigAlignment" || \
{ echo 'Problem aligning' "$ContigFile"'. Quitting.' >&2 ; exit 1 ; }

# Swap the contigs from after the references to before them, for easier visual
# inspection.
NumRefs=$(grep -e '^>' "$RefAlignment" | wc -l)
ContigsStartLine=$(awk '/^>/ {N++; if (N=='$((NumRefs+1))') {print NR; exit}}' \
"$TempRawContigAlignment")
tail -n +"$ContigsStartLine" "$TempRawContigAlignment" > "$RawContigAlignment"
head -n "$((ContigsStartLine-1))" "$TempRawContigAlignment" >> \
"$RawContigAlignment"

# Run the contig cutting & flipping code.
if Rscript "$Code_ContigCutter" "$BlastFile" "$ContigFile" "$CutContigFile"; \
then

  # The contig cutting & flipping code worked. Check if it did anything.
  "$Code_FastaFileComparer" "$RawContigFile" "$CutContigFile"
  ComparisonExitStatus=$?

  # If the contig cutting & flipping code left the contigs untouched do nothing.
  if [ $ComparisonExitStatus -eq 0 ]; then
    echo 'No cutting or flipping was necessary for the contigs.'
    # TODO: rm will no longer be necessary when the contig cutter is re-written.
    rm "$CutContigFile"

  # If the contig cutting & flipping code modified the contigs, align them.
  elif [ $ComparisonExitStatus -eq 111 ]; then
    mafft --quiet --add "$CutContigFile" "$RefAlignment" \
    > "$TempCutContigAlignment" || \
    { echo 'Problem aligning' "$CutContigFile"'. Quitting.' >&2 ; exit 1 ; }

    # Swap the contigs from after the references to before them.
    ContigsStartLine=$(awk '/^>/ {N++; if (N=='$((NumRefs+1))')'\
    '{print NR; exit}}' "$TempCutContigAlignment")
    tail -n +"$ContigsStartLine" "$TempCutContigAlignment" > \
    "$CutContigAlignment"
    head -n "$((ContigsStartLine-1))" "$TempCutContigAlignment" >> \
    "$CutContigAlignment"

  else
    echo 'Problem running' "$Code_FastaFileComparer"'. Quitting.' >&2
    exit 1
  fi  

else
  echo 'Error encountered running' "$Code_ContigCutter"
fi
