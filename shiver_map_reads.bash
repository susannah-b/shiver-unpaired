#!/usr/bin/env bash

set -u

################################################################################
# PRELIMINARIES

# Check for the right number of arguments. Assign them to variables.
NumArgsExpected=8
if [ "$#" -ne "$NumArgsExpected" ]; then
  echo "$#" 'arguments specified;' "$NumArgsExpected" 'expected. Quitting' >&2
  exit 1
fi
InitDir="$1"
ConfigFile="$2"
RawContigsFile="$3"
SID="$4"
ContigBlastFile="$5"
FastaFile="$6"
reads1="$7"
reads2="$8"

# Check InitDir exists. Remove a trailing slash, if present.
if [ ! -d "$InitDir" ]; then
  echo "$InitDir does not exist. Quitting." >&2
  exit 1
fi
InitDir=$(cd "$InitDir"; pwd)

RefList="$InitDir"/'ExistingRefNamesSorted.txt'
ExistingRefAlignment="$InitDir"/'ExistingRefAlignment.fasta'
adapters="$InitDir"/'adapters.fasta'
primers="$InitDir"/'primers.fasta'

# Source required code & check files exist
ThisDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source "$ThisDir"/'shiver_funcs.bash'
CheckFilesExist "$ConfigFile" "$reads1" "$reads2" "$RawContigsFile" \
"$ContigBlastFile" "$FastaFile" "$RefList" "$ExistingRefAlignment" "$adapters" \
"$primers"
source "$ConfigFile"

# Some files we'll create
TheRef="$SID$OutputRefSuffix"
consensus="$SID"'_MinCov_'"$MinCov1"'_'"$MinCov2.fasta"
ConsensusForGlobalAln="$SID"'_MinCov_'"$MinCov1"'_'"$MinCov2$GlobalAlnSuffix"
MappedContaminantReads="$SID$MappedContaminantReadsSuffix"
cleaned1reads="$SID$CleanedReads1Suffix"
cleaned2reads="$SID$CleanedReads2Suffix"
CoordsDict="$SID$CoordsDictSuffix"
BaseFreqs="$SID$BaseFreqsSuffix"
InsertSizeCounts="$SID$InsertSizeCountsSuffix"
################################################################################



################################################################################
# CONSTRUCT A REFERENCE, OR USE THE ONE SUPPLIED

# FastaFile should be either a single seq, which we map to as is, or else an 
# alignment of contigs to real refs.
RefIsInAlignment=true
NumSeqsInFastaFile=$(grep -e '^>' "$FastaFile" | wc -l)
if [ "$NumSeqsInFastaFile" -eq 0 ]; then
  echo 'Error: there are no sequences in' "$FastaFile". 'Quitting.' >&2
  exit 1

elif [ "$NumSeqsInFastaFile" -eq 1 ]; then

  # Try to find the sequence in FastaFile in ExistingRefAlignment.
  RefName=$(awk '/^>/ {print substr($1,2)}' "$FastaFile")
  "$Code_FindSeqsInFasta" "$ExistingRefAlignment" -g "$RefName" > \
  "$RefFromAlignment" || \
  { echo 'Could not find seq' "$RefName" 'in' "$ExistingRefAlignment"'. After' \
  'mapping we will not be able to produce a version of the consensus seq' \
  'suitable for a global alignment. Continuing.' ; RefIsInAlignment=false ; }

  # Compare the sequence in FastaFile to the one in ExistingRefAlignment.
  "$FastaFileComparer" "$RefFromAlignment" "$FastaFile"
  ComparisonExitStatus=$?
  if [ $ComparisonExitStatus -eq 111 ]; then
    echo 'Seq' "$RefName" 'differs between' "$FastaFile" 'and' \
    "$ExistingRefAlignment"'. After mapping we will not be able to produce a' \
    'version of the consensus seq suitable for a global alignment. Continuing.'
    RefIsInAlignment=false
  elif [ $ComparisonExitStatus -neq 0 ]; then
    echo 'Problem running' "$FastaFileComparer"'. Quitting.' >&2
    exit 1
  fi 

  # Set the flag appropriate for use of a real ref, when it comes to
  # coordinate translation for a global alignment.
  GlobalAlignExcisionFlag='-d'

  cp "$FastaFile" "$TheRef"
  cp "$ExistingRefAlignment" "$TempRefAlignment"

else

  ContigToRefAlignment="$FastaFile"

  #HIVcontigNames=$(awk '/^>/ {print substr($1,2)}' "$ContigToRefAlignment" | \
  #awk '/'"$SID"'/ {print}')
  #HIVcontigNames=$(awk '/^>/' "$ContigToRefAlignment" | \
  #awk '/'"$SID"'/ {printf substr($1,2) " "}')
  #HIVcontigNames=$(awk '/^>CPZ.US.85.US_Marilyn.AF103818$/ {FoundLastRef=1; next} FoundLastRef && /^>/ {printf substr($0,2) " "}' "$ContigToRefAlignment")
  #HIVcontigNames=$(awk -F, '{print $1}' "$ContigBlastFile" | sort | uniq)

  # ContigToRefAlignment should contain the same set of sequences in the input
  # existing reference alignment, plus the contigs.
  awk '/^>/ {print substr($1,2)}' "$ContigToRefAlignment" | sort > \
  "$AllSeqsInAln"
  MissingRefs=$(comm -1 -3 "$AllSeqsInAln" "$RefList")
  NumMissingRefs=$(echo $MissingRefs | wc -w)
  if [ $NumMissingRefs -gt 0 ]; then
    echo "Error: the following references from $ExistingRefAlignment are"\
    "missing from $ContigToRefAlignment: $MissingRefs. Quitting." >&2
    exit 1
  fi
  HIVcontigNames=$(comm -2 -3 "$AllSeqsInAln" "$RefList")
  NumHIVContigs=$(echo "$HIVcontigNames" | wc -w)
  if [ $NumHIVContigs -eq 0 ]; then
    echo "Error: no contigs found in $ContigToRefAlignment. Quitting" >&2
    exit 1
  fi

  # Extract just the existing references (i.e. everything but the contigs) from
  # ContigToRefAlignment, and check that they are the same as in
  # ExistingRefAlignment.
  "$Code_FindSeqsInFasta" "$ContigToRefAlignment" $HIVcontigNames -v > \
  "$TempRefAlignment" || { echo 'Problem extracting just the existing refs'\
  "from $ContigToRefAlignment. Quitting." >&2 ; exit 1 ; }
  "$Code_RemoveBlankCols" "$TempRefAlignment" > "$AlignmentForTesting" || \
  { echo "Problem removing pure-gap columns from $TempRefAlignment (which was"\
  "created by removing the contigs from $ContigToRefAlignment - that's"\
  "probably the problematic file). Quitting." >&2 ; exit 1; }
  "$Code_FastaFileComparer" "$AlignmentForTesting" "$ExistingRefAlignment"
  ComparisonExitStatus=$?
  if [ $ComparisonExitStatus -eq 111 ]; then
    echo "The reference sequences in $ContigToRefAlignment are different from"\
    "those in $ExistingRefAlignment. When modifying $ContigToRefAlignment you"\
    "should only have modified the contigs. Quitting." >&2  
    exit 1
  elif [ $ComparisonExitStatus -ne 0 ]; then
    echo 'Problem running' "$Code_FastaFileComparer"'. Quitting.' >&2
    exit 1
  fi

  # Construct the tailored ref
  "$Code_ConstructRef" "$ContigToRefAlignment" $HIVcontigNames \
  > "$GappyRefWithExtraSeq" || \
  { echo 'Failed to construct a ref from the alignment. Quitting.' >&2 ; \
  exit 1 ; }

  # Extract just the constructed ref (the first sequence)
  awk '/^>/{if(N)exit;++N;} {print;}' "$GappyRefWithExtraSeq" > "$RefWithGaps"

  # Remove any gaps from the reference
  "$Code_UngapFasta" "$RefWithGaps" > "$TheRef" || \
  { echo 'Gap stripping code failed. Quitting.' >&2 ; exit 1 ; }

  RefName=$(awk '/^>/ {print substr($1,2)}' "$TheRef")

  # Set the flag appropriate for use of a constructed ref, when it comes to
  # coordinate translation for a global alignment.
  GlobalAlignExcisionFlag='-e'

  # Create a version of the alignment of contigs to real refs, with the contigs 
  # replaced by the constructed ref, ready for coordinate translation later.
  cat "$RefWithGaps" >> "$TempRefAlignment"

fi

# Index the ref
"$smalt" index $smaltIndexOptions "$smaltIndex" "$TheRef" || \
{ echo 'Problem indexing the refererence with smalt. Quitting.' >&2 ; exit 1 ; }  


################################################################################

################################################################################
# TRIM & CLEAN READS

# Copy the reads to the working directory. Unzip them if they end in .gz.
cp "$reads1" "$reads2" .
reads1=$(basename "$reads1")
reads2=$(basename "$reads2")
if [[ "$reads1" == *.gz ]]; then
  gunzip -f "$reads1"
  reads1="${reads1%.gz}"
fi
if [[ "$reads2" == *.gz ]]; then
  gunzip -f "$reads2"
  reads2="${reads2%.gz}"
fi

# Trim reads for adapters and low-quality bases
echo 'Now trimming reads  - typically a slow step.'
java -jar "$trimmomatic" PE -threads $NumThreadsTrimmomatic \
"$reads1" "$reads2" "$reads1trim1" "$reads1trimmings" "$reads2trim1" \
"$reads2trimmings" ILLUMINACLIP:"$adapters":"$IlluminaClipParams" \
$BaseQualityParams || { echo 'Problem running trimmomatic. Quitting.' >&2 ; \
exit 1 ; }

# Trim reads for primers
$FastaqSequenceTrim "$reads1trim1" "$reads2trim1" "$reads1trim2" \
"$reads2trim2" "$primers" || \
{ echo 'Problem running fastaq. Quitting.' >&2 ; exit 1 ; }

# We'll only work with the trimmed reads now, so rename for brevity:
reads1="$reads1trim2"
reads2="$reads2trim2"

# List all the contigs and the HIV ones.
# TODO: later on we assume the blast file first field has no whitespace in it...
awk '/^>/ {print substr($1,2)}' "$RawContigsFile" | sort > "$AllContigsList"
awk -F, '{print $1}' "$ContigBlastFile" | sort | uniq > "$HIVContigsList"

# Check there are some contigs
NumContigs=$(wc -l "$AllContigsList" | awk '{print $1}')
if [ "$NumContigs" -eq 0 ]; then
  echo 'Error: there are no contigs in' "$RawContigsFile"
  echo 'Quitting.' >&2
  exit 1
fi

# Check that there aren't any contigs appearing in the blast file & missing from
# the file of contigs.
NumUnknownContigsInBlastHits=$(comm -1 -3 "$AllContigsList" "$HIVContigsList" \
| wc -l | awk '{print $1}')
if [ "$NumUnknownContigsInBlastHits" -ne 0 ]; then
  echo 'Error: the following contigs are named in' "$ContigBlastFile"\
  'but are not in' "$RawContigsFile"':'
  comm -1 -3 "$AllContigsList" "$HIVContigsList"
  echo 'Quitting.' >&2
  exit 1
fi

# Find the contaminant contigs.
ContaminantContigNames=$(comm -3 "$AllContigsList" "$HIVContigsList")
NumContaminantContigs=$(echo $ContaminantContigNames | wc -w)

# If there are no contaminant contigs, we don't need to clean.
# We create a blank mapping file to more easily keep track of the fact that 
# there are no contaminant reads in this case.
if [ "$NumContaminantContigs" -eq 0 ]; then
  echo 'There are no contaminant contigs: read cleaning unnecessary.'
  echo -n > "$MappedContaminantReads"
  mv "$reads1" "$cleaned1reads"
  mv "$reads2" "$cleaned2reads"

# We enter this scope if there are some contaminant contigs:
else

  # Make a blast database out of the contaminant contigs and the ref.
  "$Code_FindSeqsInFasta" "$RawContigsFile" $ContaminantContigNames > \
  "$RefAndContaminantContigs"
  cat "$TheRef" >> "$RefAndContaminantContigs"
  "$BlastDBcommand" -dbtype nucl -in "$RefAndContaminantContigs" \
  -input_type fasta -out "$BlastDB" || \
  { echo 'Problem creating a blast database. Quitting.' >&2 ; exit 1 ; }

  # Convert fastq to fasta.
  sed -n '1~4s/^@/>/p;2~4p' "$reads1" > "$reads1asFasta" &&
  sed -n '1~4s/^@/>/p;2~4p' "$reads2" > "$reads2asFasta" || \
  { echo 'Problem converting the reads from fastq to fasta. Quitting.' >&2 ; \
  exit 1 ; }

  # Blast the reads.
  echo 'Now blasting the reads - typically a slow step.'
  blastn -query "$reads1asFasta" -db "$BlastDB" -out "$reads1blast1" \
  -max_target_seqs 1 -outfmt \
  '10 qacc sacc sseqid evalue pident qstart qend sstart send' &&
  blastn -query "$reads2asFasta" -db "$BlastDB" -out "$reads2blast1" \
  -max_target_seqs 1 -outfmt \
  '10 qacc sacc sseqid evalue pident qstart qend sstart send' || \
  { echo 'Problem blasting' "$ContigFile"'. Quitting.' >&2 ; exit 1 ; }

  # For multiple blast hits, keep the one with the highest evalue
  # TODO: test what blast does with fasta headers that have comments in them -
  # does it include them too?
  sort -t, -k1,1 -k4,4g "$reads1blast1" | sort -t, -k1,1 -u --merge > \
  "$reads1blast2"
  sort -t, -k1,1 -k4,4g "$reads2blast1" | sort -t, -k1,1 -u --merge > \
  "$reads2blast2"

  # Find the read pairs that blast best to something other than the reference.
  "$Code_FindContaminantReadPairs" "$reads1blast2" "$reads2blast2" "$RefName" \
  "$BadReadsBaseName" && ls "$BadReadsBaseName"_1.txt \
  "$BadReadsBaseName"_2.txt > /dev/null 2>&1 || \
  { echo 'Problem finding contaminant read pairs using' \
  "$Code_FindContaminantReadPairs. Quitting." >&2 ; exit 1 ; }

  # If none of the read pairs blast better to contaminant contigs than the
  # reference, we just duplicate the original short read files.
  NumContaminantReadPairs=$(wc -l "$BadReadsBaseName"_1.txt | awk '{print $1}')
  if [ "$NumContaminantReadPairs" -eq 0 ]; then
    echo 'There are no contaminant read pairs.'
    echo -n > "$MappedContaminantReads"
    mv "$reads1" "$cleaned1reads"
    mv "$reads2" "$cleaned2reads"

  # We enter this scope if there are some read pairs that blast better to 
  # contaminant contigs than the reference.
  else

    # Sort the raw reads by name. Check every read has a mate.
    # TODO: this breaks if there's a tab in the fastq header line: grep or awk
    # for that.
    # TODO: this assumes that trimming two characters off the end fastq id field
    # gives something that matches between forward & backward reads.
    # Verify this: check that printing only the last two characters, we get only
    # one thing from each file. Move the 'unpaired' check right to the beginning.
    cat "$reads1" | paste - - - - | sort -k1,1 -t$'\t' | tr "\t" "\n" > \
    "$reads1sorted"
    cat "$reads2" | paste - - - - | sort -k1,1 -t$'\t' | tr "\t" "\n" > \
    "$reads2sorted"
    if ! cmp <(awk '{if ((NR-1)%4==0) print substr($1,2,length($1)-3)}' \
    "$reads1sorted" | sort) \
    <(awk '{if ((NR-1)%4==0) print substr($1,2,length($1)-3)}' \
    "$reads2sorted" | sort); then
      echo 'At least one read in' "$reads1" 'or' "$reads2" 'is unpaired.' \
      'Quitting.' >&2 ; exit 1 ;
    fi

    # Extract the non-contaminant read pairs
    mv "$BadReadsBaseName"_1.txt "$BadReadsBaseName"_1_unsorted.txt
    mv "$BadReadsBaseName"_2.txt "$BadReadsBaseName"_2_unsorted.txt
    sort "$BadReadsBaseName"_1_unsorted.txt > "$BadReadsBaseName"_1.txt
    sort "$BadReadsBaseName"_2_unsorted.txt > "$BadReadsBaseName"_2.txt
    "$Code_FindReadsInFastq" -v "$reads1sorted" "$BadReadsBaseName"_1.txt > \
    "$cleaned1reads" &&
    "$Code_FindReadsInFastq" -v "$reads2sorted" "$BadReadsBaseName"_2.txt > \
    "$cleaned2reads" || \
    { echo 'Problem extracting the non-contaminant reads using' \
    "$Code_FindReadsInFastq"'. Quitting.' >&2 ; exit 1 ; }

    # Map the contaminant reads to the reference, to measure how useful the
    # cleaning procedure was.
    "$Code_FindReadsInFastq" "$reads1sorted" "$BadReadsBaseName"_1.txt > \
    "$BadReadsBaseName"_1.fastq &&
    "$Code_FindReadsInFastq" "$reads2sorted" "$BadReadsBaseName"_2.txt > \
    "$BadReadsBaseName"_2.fastq || \
    { echo 'Problem extracting the contaminant reads using' \
    "$Code_FindReadsInFastq. Quitting." >&2 ; exit 1 ; }
    "$samtools" faidx "$TheRef" &&
    "$smalt" map $smaltMapOptions -o "$AllMappedContaminantReads" \
    "$smaltIndex" "$BadReadsBaseName"_1.fastq "$BadReadsBaseName"_2.fastq &&
    "$samtools" view -bS -F 4 -t "$TheRef".fai -o "$MappedContaminantReads" \
    "$AllMappedContaminantReads" || \
    { echo "Problem mapping the contaminant reads to $RefName using smalt." \
    'Quitting.' >&2 ; exit 1 ; }

  fi
fi


################################################################################
# MAP

# Do the mapping!
echo 'Now mapping - typically a slow step.'
"$smalt" map $smaltMapOptions -o "$MapOutAsSam" "$smaltIndex" "$cleaned1reads" \
"$cleaned2reads" || \
{ echo 'Smalt mapping failed. Quitting.' >&2 ; exit 1 ; }

# Convert that sam file into a bam file. Thanks Nick Croucher!
"$samtools" view -bS $samtoolsReadFlags -t "$TheRef".fai -o \
"$MapOutConversion1" "$MapOutAsSam" &&
"$samtools" sort -n "$MapOutConversion1" "$MapOutConversion2" &&
"$samtools" fixmate "$MapOutConversion2".bam "$MapOutConversion3" &&
"$samtools" sort "$MapOutConversion3" "$SID" &&
"$samtools" index "$SID.bam" || \
{ echo 'Failed to convert from sam to bam format. Quitting.' >&2 ; exit 1 ; }

# Calculate the normalised insert size distribution.
"$samtools" view "$SID.bam" | awk '{if ($9 > 0) print $9}' > "$InsertSizes1"
sort -n "$InsertSizes1" | uniq -c > "$InsertSizes2"
InsertCount=$(awk '{sum+=$1} END {print sum}' "$InsertSizes1")
awk '{print $2 "," $1 "," $1/'$InsertCount'}' "$InsertSizes2" > \
"$InsertSizeCounts"

# Generate pileup
echo 'Now calculating pileup - typically a slow step.'
"$samtools" mpileup $mpileupOptions -f "$TheRef" "$SID.bam" > "$PileupFile" || \
{ echo 'Failed to generate pileup. Quitting.' >&2 ; exit 1 ; }

# Generate base frequencies and consensuses
"$Code_AnalysePileup" "$PileupFile" "$TheRef" > "$BaseFreqs" && \
"$Code_CallConsensus" "$BaseFreqs" "$MinCov1" "$MinCov2" > \
"$consensus" || \
{ echo 'Problem analysing the pileup or calling the consensus.' >&2 ; exit 1 ; }

# Add gaps and excise unique insertions, to allow this consensus to be added to
# a global alignment with others.
if $RefIsInAlignment; then
  "$Code_MergeAlignments" "$GlobalAlignExcisionFlag" -L "$CoordsDict" \
  "$TempRefAlignment" "$consensus" > "$ConsensusForGlobalAln"
fi

# TODO: merge CoordsDict with BaseFreqs, like this?
#  ~/Dropbox\ \(Infectious\ Disease\)/chris/SeqAnal/MergeBaseFreqsAndCoords.py "$BaseFreqs" "$CoordsDict" > "$CoordsDict"_wGlobal.csv; 

# Zip the cleaned reads
gzip -f "$cleaned1reads" "$cleaned2reads"

