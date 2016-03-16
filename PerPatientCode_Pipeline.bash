#!/usr/bin/env bash

set -u

# TODO: check contig / seq names are unique? wc -l == sort | uniq | wc -l ?

################################################################################
# PRELIMINARIES

# Check for the right number of arguments. Assign them to variables.
if [ "$#" -neq 7 ]; then
  echo "$#" 'arguments specified; exactly 7 are required. Quitting' >&2
  exit 1
fi
reads1="$1"
reads2="$2"
RawContigsFile="$3"
ConfigFile="$4"
SID="$5" # Needs to match contig names in FastaFile. Otherwise just for labelling
ContigBlastFile="$6"
FastaFile="$7"
# dispense with this: just check whether FastaFile has one, or several, seq(s).
#  FastaFileType="$7"
#  if [ "$FastaFileType" == '-alignedcontigs' ]; then
#    ConstructRef=true
#  elif [ "$FastaFileType" == '-oneref' ]; then
#    ConstructRef=false
#  else
#    echo 'The seventh argument should be either -alignedcontigs or -oneref.'\
#    'Quitting' >&2
#    exit 1

source "$ConfigFile"

# Check the reads end in _1 then the specified extension
# TODO: check this wasn't needed...
#if [[ ! "$reads1" == *'_1'"$ReadExtension" ]]; then
#  echo "$reads1" 'does not terminate with _1'"$ReadExtension"'. Quitting.' >&2
#  exit 1
#fi
# TODO: define SID to be SID=$(basename "${reads1%_1$ReadExtension}") or the user arg?


################################################################################



################################################################################
# CONSTRUCT A REFERENCE, OR USE A REAL ONE

TheRef="$SID""$OutputRefSuffix"
RefWithGaps="$SID""$OutputRefWgapsSuffix"
ConsensusBasename="$SID"'_MinCov_'"$MinCoverage1"'_'"$MinCoverage2"
TestFile='dummy.fasta'

# FastaFile should either contain a single seq, which we map to as is, or else
# be an alignment of contigs to real refs.
NumSeqsInFastaFile=$(grep -e '^>' "$FastaFile" | wc -l)
if [ "$NumSeqsInFastaFile" -eq 0 ]; then
  echo 'Error: there are no sequences in' "$FastaFile". 'Quitting.' >&2
  exit 1

elif [ "$NumSeqsInFastaFile" -eq 1 ]; then

  # Try to find the sequence in FastaFile in AlignmentOfRealRefs.
  RefName=$(awk '/^>/ {print substr($1,2)}' "$FastaFile")
  RefIsInAlignment=true
  "$Code_FindSeqsInFasta" "$AlignmentOfRealRefs" -g "$RefName" > "$TestFile" ||\
  { echo 'Could not find seq' "$RefName" 'in' "$AlignmentOfRealRefs"'. After' \
  'mapping we will not be able to produce a version of the consensus seq' \
  'suitable for a global alignment. Continuing.' ; RefIsInAlignment=false ; }

  # Compare the sequence in FastaFile to the one in AlignmentOfRealRefs.
  "$FastaFileComparer" "$TestFile" "$FastaFile"
  ComparisonExitStatus=$?
  if [ $ComparisonExitStatus -eq 111 ]; then
    echo 'Seq' "$RefName" 'differs between' "$FastaFile" 'and' \
    "$AlignmentOfRealRefs"'. After mapping we will not be able to produce a' \
    'version of the consensus seq suitable for a global alignment. Continuing.'
    RefIsInAlignment=false
  elif [ $ComparisonExitStatus -eq 0 ]; then
    echo 'Problem running' "$FastaFileComparer"'. Quitting.' >&2
    exit 1
  fi 

  # Set the flag appropriate for use of a real ref, when it comes to
  # coordinate translation for a global alignment.
  GlobalAlignmentExcisionFlag='-d'

  cp "$FastaFile" "$TheRef"
  cp "$AlignmentOfRealRefs" "$TempRefAlignment"

else

  ContigToRefAlignment="$FastaFile"
  # Find the contig names in the alignment, assumed to be any sequence names
  # containing the SID. If no contigs are found, quit.
  # TODO: test this works as intended - only maches SID in the first word
  ContigNames=$(awk '/^>/ {print substr($1,2)}' "$ContigToRefAlignment" | \
  awk '/'"$SID"'/ {print}')
  #ContigNames=$(awk '/^>/' "$ContigToRefAlignment" | \
  #awk '/'"$SID"'/ {printf substr($1,2) " "}')
  #ContigNames=$(awk '/^>CPZ.US.85.US_Marilyn.AF103818$/ {FoundLastRef=1; next} FoundLastRef && /^>/ {printf substr($0,2) " "}' "$ContigToRefAlignment")
  NumContigs=$(echo $ContigNames | wc -w)
  if [ $NumContigs -eq 0 ]; then
    echo 'No contigs found in' "$ContigToRefAlignment"'. Quitting' >&2
    exit 1
  fi

  # Construct the tailored ref
  "$Code_ConstructRef" "$ContigToRefAlignment" $ContigNames \
  > "$GappyRefWithExtraSeq" || \
  { echo 'Failed to construct a ref from the alignment. Quitting.' >&2 ; \
  exit 1 ; }

  # Extract just the constructed ref (the first sequence)
  awk '/^>/{if(N)exit;++N;} {print;}' "$GappyRefWithExtraSeq" \
  > "$RefWithGaps"

  # Remove any gaps from the reference
  "$Code_RemoveGaps" "$RefWithGaps" > "$TheRef" || \
  { echo 'Gap stripping code failed. Quitting.' >&2 ; exit 1 ; }

  RefName=$(awk '/^>/ {print substr($1,2)}' "$TheRef")

  # Set the flag appropriate for use of a constructed ref, when it comes to
  # coordinate translation for a global alignment.
  GlobalAlignmentExcisionFlag='-e'

  # In the alignment of contigs to real refs, replace the contigs by the
  # constructed ref ready for coordinate translation later.
  # Do this by first extracting all sequences except the contigs, then adding
  # the constructed ref with gaps.
  "$Code_FindSeqsInFasta" "$ContigToRefAlignment" $ContigNames -v > \
  "$TempRefAlignment" || { echo 'Problem extracting just the real refs from '\
  "$ContigToRefAlignment"'. Quitting.' >&2 ; exit 1 ; }  
  cat "$RefWithGaps" >> "$TempRefAlignment"

fi

# Index the ref
"$smalt" index $smaltIndexOptions "$smaltIndex" "$TheRef" || \
{ echo 'Problem indexing the refererence with smalt. Quitting.' >&2 ; exit 1 ; }  


################################################################################


################################################################################
# REMOVE CONTAMINANT READS

# Copy some files to the working directory.
cp "$RawContigsFile" "$ContigBlastFile" .
RawContigsFile=$(basename "$RawContigsFile")
ContigBlastFile=$(basename "$ContigBlastFile")

# List all the contigs and the HIV ones.
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

MappedContaminantReads="$SID""$MappedContaminantReadsSuffix"

# Find the contaminant contigs.
comm -3 "$AllContigsList" "$HIVContigsList" > "$ContaminantContigsList"
NumContaminantContigs=$(wc -l "$ContaminantContigsList" | awk '{print $1}')

# Copy the reads to the working directory. Unzip them if they end in .gz.
cp "$reads1" "$reads2" .
reads1=$(basename "$reads1")
reads2=$(basename "$reads2")
if [[ "$ReadExtension" == *.gz ]]; then
  gunzip -f "$reads1"
  gunzip -f "$reads2"
  reads1="${reads1%.gz}"
  reads2="${reads2%.gz}"
  ReadExtension="${ReadExtension%.gz}"
fi

# If there are no contaminant contigs, we don't need to clean.
# We create a blank mapping file to indicate that there wasn't a failure to
# map, it was just unnecessary.
if [ "$NumContaminantContigs" -eq 0 ]; then
  echo 'There are no contaminant contigs. The cleaned reads = the orginal reads'
  echo -n > "$MappedContaminantReads"
  #cp "$reads1" "$SID"'_cleaned_1'"$ReadExtension"
  #cp "$reads2" "$SID"'_cleaned_2'"$ReadExtension"
  #if [[ "$ReadExtension" == *.gz ]]; then
  #  mv "$reads1" "$cleaned1reads"
  #  mv "$reads2" "$cleaned2reads"
  #else
  #  cleaned1reads="$cleaned1reads".gz
  #  cleaned2reads="$cleaned2reads".gz
  #  gzip -f "$reads1" "$reads2" 
  #  mv "$reads1".gz "$cleaned1reads"
  #  mv "$reads2".gz "$cleaned2reads"
  #fi
  cleaned1reads="$reads1"
  cleaned2reads="$reads2"

# We enter this scope if there are some contaminant contigs:
else

  # Make a blast database out of the contaminant contigs and the ref.
  "$Code_FindSeqsInFasta" "$RawContigsFile" $(cat "$ContaminantContigsList" | xargs echo) \
  > "$RefAndContaminantContigs"
  cat "$TheRef" >> "$RefAndContaminantContigs"
  "$BlastDBcommand" -dbtype nucl -in "$RefAndContaminantContigs" \
  -input_type fasta -out "$BlastDB" || \
  { echo 'Problem creating a blast database. Quitting.' >&2 ; exit 1 ; }

  # Convert fastq to fasta.
  reads1asFasta="${reads1%$ReadExtension}".fasta
  reads2asFasta="${reads2%$ReadExtension}".fasta
  sed -n '1~4s/^@/>/p;2~4p' "$reads1" > "$reads1asFasta" &&
  sed -n '1~4s/^@/>/p;2~4p' "$reads2" > "$reads2asFasta" || \
  { echo 'Problem converting the reads from fastq to fasta. Quitting.' >&2 ; \
  exit 1 ; }

  # Blast the reads.
  reads1blast="${reads1%$ReadExtension}".blast
  reads2blast="${reads2%$ReadExtension}".blast
  blastn -query "$reads1asFasta" -db "$BlastDB" -out "$reads1blast".temp \
  -max_target_seqs 1 -outfmt \
  '10 qacc sacc sseqid evalue pident qstart qend sstart send' &&
  blastn -query "$reads2asFasta" -db "$BlastDB" -out "$reads2blast".temp \
  -max_target_seqs 1 -outfmt \
  '10 qacc sacc sseqid evalue pident qstart qend sstart send' || \
  { echo 'Problem blasting' "$ContigFile"'. Quitting.' >&2 ; exit 1 ; }

  # For multiple blast hits, keep the one with the highest evalue
  sort -t, -k1,1 -k4,4g "$reads1blast".temp | sort -t, -k1,1 -u --merge > \
  "$reads1blast"
  sort -t, -k1,1 -k4,4g "$reads2blast".temp | sort -t, -k1,1 -u --merge > \
  "$reads2blast"

  # Find the read pairs that blast best to something other than the reference.
  "$Code_FindContaminantReadPairs" "$reads1blast" "$reads2blast" "$RefName" \
  "$BadReadsBaseName" && ls "$BadReadsBaseName"_1.txt \
  "$BadReadsBaseName"_2.txt > /dev/null 2>&1 || \
  { echo 'Problem finding contaminant read pairs using' \
  "$Code_FindContaminantReadPairs"'. Quitting.' >&2 ; exit 1 ; }

  # If none of the read pairs blast better to contaminant contigs than the
  # reference, we just duplicate the original short read files.
  NumContaminantReadPairs=$(wc -l "$BadReadsBaseName"_1.txt | awk '{print $1}')
  if [ "$NumContaminantReadPairs" -eq 0 ]; then
    echo 'There are no contaminant read pairs. Duplicating the short read files.'
    echo -n > "$MappedContaminantReads"
    #cleaned1reads="$cleaned1reads".gz
    #cleaned2reads="$cleaned2reads".gz
    #gzip -f "$reads1" "$reads2" 
    #mv "$reads1".gz "$cleaned1reads"
    #mv "$reads2".gz "$cleaned2reads"
    cleaned1reads="$reads1"
    cleaned2reads="$reads2"

  # We enter this scope if there are some read pairs that blast better to 
  # contaminant contigs than the reference.
  else

    cleaned1reads="$SID"'_cleaned_1'"$ReadExtension"
    cleaned2reads="$SID"'_cleaned_2'"$ReadExtension"

    # Sort the raw reads by name. Check every read has a mate.
    # TODO: this breaks if there's a tab in the fastq header line: grep or awk
    # for that.
    # TODO: this assumes that trimming two characters off the end fastq id field
    # gives something that matches between forward & backward reads.
    # Verify this: check that printing only the last two characters, we get only
    # one thing from each file. Move the 'unpaired' check right to the beginning.
    cat "$reads1" | paste - - - - | sort -k1,1 -t$'\t' | tr "\t" "\n" > \
    "$reads1".sorted
    cat "$reads2" | paste - - - - | sort -k1,1 -t$'\t' | tr "\t" "\n" > \
    "$reads2".sorted
    if ! cmp <(awk '{if ((NR-1)%4==0) print substr($1,2,length($1)-3)}' \
    "$reads1".sorted | sort) \
    <(awk '{if ((NR-1)%4==0) print substr($1,2,length($1)-3)}' \
    "$reads2".sorted | sort); then
      echo 'At least one read in' "$reads1" 'or' "$reads2" 'is unpaired.' \
      'Quitting.' >&2 ; exit 1 ;
    fi

    # Extract the non-contaminant read pairs
    mv "$BadReadsBaseName"_1.txt "$BadReadsBaseName"_1_unsorted.txt
    mv "$BadReadsBaseName"_2.txt "$BadReadsBaseName"_2_unsorted.txt
    sort "$BadReadsBaseName"_1_unsorted.txt > "$BadReadsBaseName"_1.txt
    sort "$BadReadsBaseName"_2_unsorted.txt > "$BadReadsBaseName"_2.txt
    "$Code_FindReadsInFastq" -v "$reads1".sorted "$BadReadsBaseName"_1.txt > \
    "$cleaned1reads" &&
    "$Code_FindReadsInFastq" -v "$reads2".sorted "$BadReadsBaseName"_2.txt > \
    "$cleaned2reads" || \
    { echo 'Problem extracting the non-contaminant reads using' \
    "$Code_FindReadsInFastq"'. Quitting.' >&2 ; exit 1 ; }

    # Zip the cleaned reads
    #gzip -f "$cleaned1reads" "$cleaned2reads"
    #cleaned1reads="$cleaned1reads".gz
    #cleaned2reads="$cleaned2reads".gz

    # Map the contaminant reads to the reference, to measure how useful the
    # cleaning procedure was.
    "$Code_FindReadsInFastq" "$reads1".sorted "$BadReadsBaseName"_1.txt > \
    "$BadReadsBaseName"_1.fastq &&
    "$Code_FindReadsInFastq" "$reads2".sorted "$BadReadsBaseName"_2.txt > \
    "$BadReadsBaseName"_2.fastq || \
    { echo 'Problem extracting the contaminant reads using' \
    "$Code_FindReadsInFastq"'. Quitting.' >&2 ; exit 1 ; }
    "$samtools" faidx "$TheRef" &&
    "$smalt" map $smaltMapOptions -o "$AllMappedContaminantReads" \
    "$smaltIndex" "$BadReadsBaseName"_1.fastq "$BadReadsBaseName"_2.fastq &&
    "$samtools" view -bS -F 4 -t "$TheRef".fai -o "$MappedContaminantReads" \
    "$AllMappedContaminantReads" || \
    { echo 'Problem mapping the contaminant reads to' "$RefName" 'using' \
    'smalt. Quitting.' >&2 ; exit 1 ; }

  fi
fi


################################################################################


# Do the mapping!
"$smalt" map $smaltMapOptions -o "$SID".sam "$smaltIndex" "$cleaned1reads" \
"$cleaned2reads" || \
{ echo 'Smalt mapping failed. Quitting.' >&2 ; exit 1 ; }

# Convert that sam file into a bam file.
"$samtools" view -bS $samtoolsReadFlags -t "$TheRef".fai -o "$SID"_step1.bam \
"$SID".sam &&
"$samtools" sort -n "$SID"_step1.bam "$SID"_step2 &&
"$samtools" fixmate "$SID"_step2.bam "$SID"_step3.bam &&
"$samtools" sort "$SID"_step3.bam "$SID" &&
"$samtools" index "$SID""$OutputBamSuffix" || \
{ echo 'Failed to convert from sam to bam format. Quitting.' >&2 ; exit 1 ; }

"$samtools" view "$SID""$OutputBamSuffix" | awk '{if ($9 > 0) print $9}' > \
"$SID"_InsertSizes.txt
sort -n "$SID"_InsertSizes.txt | uniq -c > "$SID"_InsertSizeCounts_temp.txt
InsertCount=$(awk '{sum+=$1} END {print sum}' "$SID"_InsertSizes.txt)

awk '{print $2 "," $1 "," $1/'$InsertCount'}' \
"$SID"_InsertSizeCounts_temp.txt > "$SID""$InsertSizeCountsSuffix"

# Generate pileup
"$samtools" mpileup $mpileupOptions -f "$TheRef" "$SID""$OutputBamSuffix" > \
"$SID".pileup || { echo 'Failed to generate pileup. Quitting.' >&2 ; exit 1 ; }

# Generate base frequencies and consensuses
"$Code_AnalysePileup" "$SID".pileup "$TheRef" > "$SID""$BaseFreqsSuffix" && \
"$Code_CallConsensus" "$SID""$BaseFreqsSuffix" "$MinCoverage1" "$MinCoverage2" > \
"$ConsensusBasename""$FastaExtension" || \
{ echo 'Problem analysing the pileup or calling the consensus.' >&2 ; exit 1 ; }

# Add gaps and excise unique insertions, to allow this consensus to be added to
# a global alignment with others.
if $RefIsInAlignment; then
  "$Code_MergeAlignments" "$GlobalAlignmentExcisionFlag" -L \
  "$SID""$CoordsDictSuffix" "$TempRefAlignment" \
  "$ConsensusBasename""$FastaExtension" > \
  "$ConsensusBasename"'_ForGlobalAln'"$FastaExtension"
fi


