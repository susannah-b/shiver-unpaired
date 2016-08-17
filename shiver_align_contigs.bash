#!/usr/bin/env bash

UsageInstructions=$(echo '
Arguments for this script:
(1) the initialisation directory you created using the shiver_init.bash command;
(2) the configuration file, containing all your parameter choices etc.;
(3) a fasta file of contigs (output from processing the short reads with an
assembly program);
(4) A sample ID ("SID") used for naming the output from this script (a sensible
choice might be the contig file name minus its path and extension).
If this script completes successfully it will produce a .blast file (detailing
the contigs that have blast hits). If the .blast file is not empty it will
produce a fasta file of those contigs with hits and another fasta file of
these contigs aligned to the refereces; if cutting and/or reversing of contigs
is necessary, two more fasta files are produced - the cut/reversed contigs on
their own and also aligned to references (i.e. there will be two files of the
contigs on their own and two files of the contigs aligned to references).
')

set -u

# Check for the right number of arguments. Assign them to variables.
NumArgsExpected=4
if [ "$#" -ne "$NumArgsExpected" ]; then
  echo $UsageInstructions
  echo "$#" 'arguments specified;' "$NumArgsExpected" 'expected. Quitting' >&2
  exit 1
fi
InitDir="$1"
ConfigFile="$2"
ContigFile="$3"
SID="$4"

# Check InitDir exists. Remove a trailing slash, if present.
if [ ! -d "$InitDir" ]; then
  echo "$InitDir does not exist. Quitting." >&2
  exit 1
fi
InitDir=$(cd "$InitDir"; pwd)

BlastDatabase="$InitDir/ExistingRefsBlastDatabase"
RefAlignment="$InitDir/ExistingRefAlignment.fasta"

# Source required code & check files exist
ThisDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source "$ThisDir"/'shiver_funcs.bash'
CheckFilesExist "$ContigFile" "$RefAlignment"
source "$ConfigFile"

# Check that there are some contigs, that their IDs are unique, and that their
# IDs don't contain commas.
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
if [[ "$ContigNames" == *","* ]]; then
  echo "Contig names must not contain commas. Quitting."
  exit 1
fi

# The names for output files we'll produce.
BlastFile="$SID"'.blast'
RawContigAlignment="$SID"'_raw_wRefs.fasta'
CutContigAlignment="$SID"'_cut_wRefs.fasta'

# Blast the contigs
blastn -query "$ContigFile" -db "$BlastDatabase" -out "$BlastFile" \
-max_target_seqs 1 -outfmt \
'10 qseqid sseqid evalue pident qlen qstart qend sstart send' || \
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
"$Code_FindSeqsInFasta" "$ContigFile" $HIVcontigNames > "$RawContigFile1" || \
{ echo 'Problem extracting the HIV contigs. Quitting.' >&2 ; exit 1 ; }

# ... and align them. 
OldMafft=false
SwapContigsToTop=true
AlignContigsToRefs "$mafft" '--quiet' "$RawContigFile1" "$RefAlignment" \
"$RawContigAlignment" "$SwapContigsToTop" "$OldMafft"

# Run the contig cutting & flipping code
"$Code_CorrectContigs" "$BlastFile" -C "$ContigFile" -O "$CutContigFile" || \
{ echo "Problem encountered running $Code_CorrectContigs. Quitting." ; exit 1; }

# If the contigs needed cutting and/or flipping: align the modified contigs.
if [[ -f "$CutContigFile" ]]; then
  AlignContigsToRefs "$mafft" '--quiet' "$CutContigFile" "$RefAlignment" \
  "$CutContigAlignment" "$SwapContigsToTop" "$OldMafft"
fi

