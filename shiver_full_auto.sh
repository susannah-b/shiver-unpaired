#!/usr/bin/env bash

set -u
set -o pipefail

UsageInstructions=$(echo '
Arguments for this script:
(1) the initialisation directory you created using the shiver_init.sh command;
(2) the configuration file, containing all your parameter choices etc.;
(3) the forward reads;
(4) the reverse reads;
(5) a fasta file of contigs (output from processing the short reads with an
assembly program),
(6) A sample ID ("SID") used for naming the output from this script (a sensible
choice might be the contig or reads file name minus the path and extension).
')

# Check for the right number of arguments. Assign them to variables.
NumArgsExpected=6
if [ "$#" -ne "$NumArgsExpected" ]; then
  echo $UsageInstructions
  echo "$#" 'arguments specified;' "$NumArgsExpected" 'expected. Quitting.' >&2
  exit 1
fi
InitDir="$1"
ConfigFile="$2"
reads1="$3"
reads2="$4"
ContigFile="$5"
SID="$6"

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
source "$ThisDir"/'shiver_funcs.sh'
CheckFilesExist "$ContigFile" "$RefAlignment" "$reads1" "$reads2"
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
  echo "Contig names must not contain commas. Quitting." >&2
  exit 1
fi

# The names for output files we'll produce.
BlastFile="$SID.blast"
BestContigToRefAlignment="$SID$BestContigToRefAlignmentSuffix"
TheRef="$SID$OutputRefSuffix"

# Blast the contigs
"$BlastNcommand" -query "$ContigFile" -db "$BlastDatabase" -out "$BlastFile" \
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

# Extract those contigs that have a blast hit.
HIVcontigNames=$(awk -F, '{print $1}' "$BlastFile" | sort | uniq)
"$Code_FindSeqsInFasta" "$ContigFile" $HIVcontigNames > "$RawContigFile1" || \
{ echo 'Problem extracting the HIV contigs. Quitting.' >&2 ; exit 1 ; }

# Run the contig cutting & flipping code. Check it works, and quit if it thinks
# the contigs need correcting.
"$Code_CorrectContigs" "$BlastFile" --min-hit-frac "$MinContigHitFrac" || \
{ echo "The contigs in $ContigFile appear to need correcting. (Or a problem" \
"was encountered running $Code_CorrectContigs.) Fully automatic processing" \
"not possible. Quitting." >&2 ; exit 1 ; }

# Prepare to iterate through all existing references
mkdir "$ContigAlignmentsToRefsDir" || \
{ echo "Problem making the directory" "$ContigAlignmentsToRefsDir. (NB it" \
"should not exist already.) Quitting." >&2; exit 1; }
echo -n '' > "$RefMatchLog"
NumRefs=$(ls "$InitDir"/'IndividualRefs'/*'.fasta' | wc -l)
CurrentRef=0
OldMafft=false

# Iterate through all existing references
for ref in "$InitDir"/'IndividualRefs'/*'.fasta'; do

  CurrentRef=$((CurrentRef+1))
  echo "Aligning contigs to ref $CurrentRef of $NumRefs."

  # Make one alignment with regular mafft...
  cat "$ref" "$RawContigFile1" > "$ContigsWith1ref"
  AlnFile1="$ContigAlignmentsToRefsDir"/'temp_1_'$(basename "$ref")
  "$mafft" --quiet "$ContigsWith1ref" > "$AlnFile1" || \
  { echo 'Problem aligning' "$ContigsWith1ref"'. Quitting.' >&2 ; exit 1 ; }
  NumWordsInAlnFileName=$(echo $AlnFile1 | wc -w)
  if [[ $NumWordsInAlnFileName -gt 1 ]]; then
    echo "Error: unexpected whitespace in file $AlnFile1. Quitting." >&2
    exit 1
  fi
  RefSimilarityScore1=$("$Code_ConstructRef" -P "$AlnFile1" $HIVcontigNames | \
  awk '{print $1}')
  echo $RefSimilarityScore1 "$AlnFile1" >> "$RefMatchLog"

  # ...and try to make a second with mafft --addfragments
  if $OldMafft; then
    continue
  fi
  AlnFile2="$ContigAlignmentsToRefsDir"/'temp_2_'$(basename "$ref")
  "$mafft" --quiet --addfragments "$RawContigFile1" "$ref" > "$AlnFile2" || \
  { echo "Warning: it looks like you're running an old version of mafft: the" \
  "--addfragments option doesn't work. That option can be very helpful for" \
  "correctly aligning contigs, and we advise you to update your mafft." \
  "Continuing without using that option."
  OldMafft=true
  rm "$AlnFile2"
  continue ; }
  RefSimilarityScore2=$("$Code_ConstructRef" -P "$AlnFile2" $HIVcontigNames | \
  awk '{print $1}')
  echo $RefSimilarityScore2 "$AlnFile2" >> "$RefMatchLog"

  # TODO: mafft options, including gap extension parameter!
done

# Find the reference most closely matching the contigs
read BestRefScore BestRefFile <<<$(sort -nrk1,1 "$RefMatchLog" | head -1)
mv "$BestRefFile" "$BestContigToRefAlignment"

# Check that the gappiest contig is not too gappy in the best contig-to-ref
# alignment.
LargestGapFrac=$("$Code_ConstructRef" -S1 "$BestContigToRefAlignment" \
$HIVcontigNames | sort -nrk2,2 | head -1 | awk '{print $2}')
if (( $(echo "$LargestGapFrac > $MaxContigGappiness" | bc -l) )); then
  echo "The gappiest contig in the best contig-to-reference alignment," \
  "$BestContigToRefAlignment, has gap fraction $LargestGapFrac which is larger"\
  "than your stated maximum of $MaxContigGappiness. Try running the 'nearly" \
  "automatic' version of shiver on this sample. Quitting." >&2
  exit 1
fi

# Extract a gapless form of the best ref from its alignment to the contigs.
#"$Code_FindSeqsInFasta" "$BestContigToRefAlignment" $HIVcontigNames -v -g > \
#"$TheRef" || \
#{ echo "Problem extracting the ref from $BestRefFile. Quitting." >&2; exit 1 ; }
#NumSeqs=$(grep -e '^>' "$TheRef" | wc -l)
#if [[ $NumSeqs -ne 1 ]]; then
#  echo "Error: $TheRef contains $NumSeqs; it should contain exactly 1." \
#  "Quitting." >&2
#  exit 1
#fi

# Construct the tailored ref
"$Code_ConstructRef" "$BestContigToRefAlignment" $HIVcontigNames \
> "$GappyRefWithExtraSeq" || \
{ echo 'Failed to construct a ref from the alignment. Quitting.' >&2; exit 1 ; }

# Extract just the constructed ref (the first sequence)
awk '/^>/{if(N)exit;++N;} {print;}' "$GappyRefWithExtraSeq" > "$RefWithGaps"

# Remove any gaps from the reference
"$Code_UngapFasta" "$RefWithGaps" > "$TheRef" || \
{ echo 'Gap stripping code failed. Quitting.' >&2 ; exit 1 ; }

# Map!
"$ThisDir"/'shiver_map_reads.sh' "$InitDir" "$ConfigFile" "$ContigFile" \
"$SID" "$BlastFile" "$TheRef" "$reads1" "$reads2" || \
{ echo "Problem calling shiver_map_reads.sh. Quitting." >&2; exit 1 ; }



