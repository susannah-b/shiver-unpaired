#!/usr/bin/env bash

set -u
set -o pipefail

ThisDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ToolsDir="$ThisDir"/tools
Code_AlignToConsensus="$ToolsDir/AlignMoreSeqsToPairWithMissingCoverage.py"
Code_AnalysePileup="$ToolsDir/AnalysePileup.py"
Code_CallConsensus="$ToolsDir/CallConsensus.py"
Code_CheckFastaFileEquality="$ToolsDir/CheckFastaFileEquality.py"
Code_ConstructRef="$ToolsDir/ConstructBestRef.py"
Code_CorrectContigs="$ToolsDir/CorrectContigs.py"
Code_FillConsensusGaps="$ToolsDir/FillConsensusGaps.py"
Code_FindContaminantReadPairs="$ToolsDir/FindContaminantReadPairs.py"
Code_FindReadsInFastq="$ToolsDir/FindNamedReadsInSortedFastq.py"
Code_FindSeqsInFasta="$ToolsDir/FindSeqsInFasta.py"
Code_MergeAlignments="$ToolsDir/MergeAlignments.py"
Code_RemoveBlankCols="$ToolsDir/RemoveBlankColumns.py"
Code_SplitFasta="$ToolsDir/SplitFasta.py"
Code_UngapFasta="$ToolsDir/UngapFasta.py"
Code_MergeBaseFreqsAndCoords="$ToolsDir/MergeBaseFreqsAndCoords.py"
Code_SplitGappyContigs="$ToolsDir/SplitGappyContigs.py"

# For quitting if files don't exist.
function CheckFilesExist {
  for argument in "$@"; do
    if [ ! -f "$argument" ]; then
      echo "$argument" 'does not exist or is not a regular file. Quitting.' >&2
      exit 1
    fi
  done
}

function AlignContigsToRefs {

  Aligner=$1
  AlignerOptions=$2
  ContigFile=$3
  ThisRefAlignment=$4
  OutputAlignedContigFile=$5
  SwapContigsToTopArg=$6
  OldMafftArg=$7
  ContigNames=$(awk '/^>/ {print substr($1,2)}' "$ContigFile")

  "$Aligner" $AlignerOptions --add "$ContigFile" "$ThisRefAlignment" > \
  "$TempContigAlignment1" || \
  { echo 'Problem aligning' "$ContigFile"'.' >&2 ; return 1 ; }
  MaxContigGappiness1=$("$Code_ConstructRef" -S1 "$TempContigAlignment1" \
  $ContigNames | sort -nrk2,2 | head -1 | awk '{print $2}') || { echo 'Problem'\
  "analysing $TempContigAlignment1 (i.e. the output from aligning $ContigFile"\
  "and $ThisRefAlignment) with $Code_ConstructRef." >&2 ; return 1 ; }
  BestContigAlignment="$TempContigAlignment1"

  # If the function was called with OldMafftArg=true, don't bother trying the
  # addfragments option. Otherwise, try. If it doesn't work, set OldMafft=true
  # ready for the next call of this function; if it does work, use the least
  # gappy alignment.
  if ! $OldMafftArg; then
    "$Aligner" $AlignerOptions --addfragments "$ContigFile" \
    "$ThisRefAlignment" > "$TempContigAlignment2" || \
    { echo "Warning: it looks like you're running an old version of mafft: the"\
    "--addfragments option doesn't work. That option can be very helpful for" \
    "correctly aligning contigs, and we advise you to update your mafft." \
    "Continuing without using that option."
    OldMafft=true ; }
    NumLinesInAln=$(wc -l "$TempContigAlignment2" | awk '{print $1}')
    if ! $OldMafft; then
      MaxContigGappiness2=$("$Code_ConstructRef" -S1 "$TempContigAlignment2" \
      $ContigNames | sort -nrk2,2 | head -1 | awk '{print $2}') || { echo \
      "Problem analysing $TempContigAlignment2 (i.e. the output from aligning "\
      "$ContigFile and $ThisRefAlignment) with $Code_ConstructRef." \
      >&2 ; return 1 ; }
      if (( $(echo "$MaxContigGappiness2 < $MaxContigGappiness1" | bc -l) )); 
      then
        BestContigAlignment="$TempContigAlignment2"
        echo 'Info: the --addfragments mafft option produced a less'\
        "gappy alignment than --add did when aligning $ContigFile. Using the"\
        "--addfragments result."
      else
        echo 'Info: the --addfragments mafft option did not produce a less'\
        "gappy alignment than --add did when aligning $ContigFile. Using the"\
        "--add result."
      fi
    fi
  fi

  # Swap the contigs from after the references to before them, for easier visual
  # inspection.
  if [[ "$SwapContigsToTopArg" == "true" ]]; then
    NumRefs=$(grep -e '^>' "$ThisRefAlignment" | wc -l)
    ContigsStartLine=$(awk '/^>/ {N++; if (N=='$((NumRefs+1))') {print NR; exit}}' \
    "$BestContigAlignment")
    tail -n +"$ContigsStartLine" "$BestContigAlignment" > \
    "$OutputAlignedContigFile"
    head -n "$((ContigsStartLine-1))" "$BestContigAlignment" >> \
    "$OutputAlignedContigFile"
  else
    mv "$BestContigAlignment" "$OutputAlignedContigFile"
  fi

}

function CheckReadNames {

  ReadFile=$1
  # The second argument is 1 for forward reads or 2 for reverse reads.
  OneOrTwo=$2

  # Check all read seq ids end in /1 or /2 as needed.
  suffix=$(awk '{if ((NR-1)%4==0) print substr($1,length($1)-1,
  length($1))}' "$ReadFile" | sort | uniq)
  if [[ "$suffix" != '/'"$OneOrTwo" ]]; then
    echo "Found at least one read in $ReadFile whose sequence ID does not end"\
    'in "\'"$OneOrTwo"'".' >&2
    return 1
  fi

  # Check none of the lines with read IDs contain tabs
  NumNameLinesWithTabs=$(awk '{if ((NR-1)%4==0 && gsub("\t","\t",$0) > 0)
  print}' "$ReadFile" | wc -l)
  if [[ $NumNameLinesWithTabs -ne 0 ]]; then
    echo "The following lines in $ReadFile contain tabs:"
    awk '{if ((NR-1)%4==0 && gsub("\t","\t",$0) > 0) print}' "$ReadFile"
    echo 'To remove contaminant reads, we require there to be no tabs in the' \
    'sequence ID lines of fastq files.' >&2
    return 1
  fi

  # Check all read names are unique
  NumDuplicatedReadNames=$(awk '{if ((NR-1)%4==0) print substr($1,1)}'\
  "$ReadFile" | sort | uniq -d | wc -l)
  if [[ $NumDuplicatedReadNames -ne 0 ]]; then
    echo "The following read names are duplicated in $ReadFile:" >&2
    awk '{if ((NR-1)%4==0) print substr($1,1)}' "$ReadFile" | sort | uniq -d >&2
    echo "Reads should be uniquely named." >&2
    return 1
  fi

}


function sam_to_bam {

  # Check for the right number of args
  ExpectedNumArgs=3
  if [[ "$#" -ne "$ExpectedNumArgs" ]]; then
    echo "sam_to_bam function called with $# args; expected $ExpectedNumArgs."\
    "Quitting." >&2
    return 1
  fi

  # Assign the args
  InSam=$1
  LocalRefFAIindex=$2
  OutBam=$3

  # Thanks to Nick Croucher for these steps.
  "$samtools" view -bS $samtoolsReadFlags -t "$LocalRefFAIindex" -o \
  "$MapOutConversion1".bam "$InSam" &&
  "$samtools" sort -n "$MapOutConversion1".bam -o "$MapOutConversion2".bam -T \
  "$SamtoolsSortFile" &&
  "$samtools" fixmate "$MapOutConversion2".bam "$MapOutConversion3".bam &&
  "$samtools" sort "$MapOutConversion3".bam -o "$OutBam" -T \
  "$SamtoolsSortFile" ||
  { echo 'Failed to convert from sam to bam format.' >&2 ; return 1 ; }
}

function map_with_smalt {

  # Check for the right number of args
  ExpectedNumArgs=4
  if [[ "$#" -ne "$ExpectedNumArgs" ]]; then
    echo "map_with_smalt function called with $# args; expected"\
    "$ExpectedNumArgs. Quitting." >&2
    return 1
  fi

  # Assign the args
  ReadsToMap1=$1
  ReadsToMap2=$2
  LocalRef=$3
  OutFileAsBam=$4

  # Index the ref
  "$smalt" index $smaltIndexOptions "$smaltIndex" "$LocalRef" ||
  { echo 'Problem indexing the refererence with smalt.' >&2 ;
  return 1 ; }

  # Do the mapping!
  echo "Now mapping using smalt with options \"$smaltMapOptions\". Typically a"\
  "slow step."
  "$smalt" map $smaltMapOptions -o "$MapOutAsSam" "$smaltIndex" \
  "$ReadsToMap1" "$ReadsToMap2" || \
  { echo 'Smalt mapping failed.' >&2 ; return 1 ; }

  # Make the reference's .fai index if needed.
  LocalRefFAIindex="$LocalRef".fai
  if [[ ! -f "$LocalRefFAIindex" ]]; then
    "$samtools" faidx "$LocalRef" && ls "$LocalRef".fai > /dev/null ||
    { echo 'Problem indexing the refererence with samtools. Quitting.' >&2 ; 
    return 1 ; }
  fi

  sam_to_bam "$MapOutAsSam" "$LocalRefFAIindex" "$OutFileAsBam" ||
  { echo 'Problem converting from sam to bam format.' >&2 ; return 1 ; }
}

function map_with_bwa_mem {

  # Check for the right number of args
  ExpectedNumArgs=4
  if [[ "$#" -ne "$ExpectedNumArgs" ]]; then
    echo "map_with_bwa_mem function called with $# args; expected"\
    "$ExpectedNumArgs. Quitting." >&2
    return 1
  fi

  # Assign the args
  ReadsToMap1=$1
  ReadsToMap2=$2
  LocalRef=$3
  OutFileAsBam=$4

  # Index the ref
  "$bwa" index "$LocalRef" ||
  { echo 'Problem indexing the refererence with bwa.' >&2 ;
  return 1 ; }

  # Do the mapping!
  echo "Now mapping using bwa mem with options \"$bwaOptions\". Typically a"\
  "slow step."
  "$bwa" mem "$LocalRef" "$ReadsToMap1" "$ReadsToMap2" $bwaOptions > \
  "$MapOutAsSam" || { echo 'bwa mem mapping failed.' >&2 ; return 1 ; }

  # Make the reference's .fai index if needed.
  LocalRefFAIindex="$LocalRef".fai
  if [[ ! -f "$LocalRefFAIindex" ]]; then
    "$samtools" faidx "$LocalRef" && ls "$LocalRef".fai > /dev/null ||
    { echo 'Problem indexing the refererence with samtools. Quitting.' >&2 ; 
    return 1 ; }
  fi

  sam_to_bam "$MapOutAsSam" "$LocalRefFAIindex" "$OutFileAsBam" ||
  { echo 'Problem converting from sam to bam format.' >&2 ; return 1 ; }
}

function map_with_bowtie {

  # Check for the right number of args
  ExpectedNumArgs=4
  if [[ "$#" -ne "$ExpectedNumArgs" ]]; then
    echo "map_with_bowtie function called with $# args; expected"\
    "$ExpectedNumArgs. Quitting." >&2
    return 1
  fi

  # Assign the args
  ReadsToMap1=$1
  ReadsToMap2=$2
  LocalRef=$3
  OutFileAsBam=$4

  # Index the ref
  "$bowtie2_build" "$LocalRef" "$bowtieIndex" ||
  { echo 'Problem indexing the refererence with bowtie2.' >&2 ;
  return 1 ; }

  # Do the mapping!
  echo "Now mapping using bowtie2 with options \"$bowtieOptions\". Typically a"\
  "slow step."
  "$bowtie2" -x "$bowtieIndex" -1 "$ReadsToMap1" -2 "$ReadsToMap2" -S \
  "$MapOutAsSam" $bowtieOptions || \
  { echo 'bowtie2 mapping failed.' >&2 ; return 1 ; }

  # Make the reference's .fai index if needed.
  LocalRefFAIindex="$LocalRef".fai
  if [[ ! -f "$LocalRefFAIindex" ]]; then
    "$samtools" faidx "$LocalRef" && ls "$LocalRef".fai > /dev/null ||
    { echo 'Problem indexing the refererence with samtools. Quitting.' >&2 ; 
    return 1 ; }
  fi

  sam_to_bam "$MapOutAsSam" "$LocalRefFAIindex" "$OutFileAsBam" ||
  { echo 'Problem converting from sam to bam format.' >&2 ; return 1 ; }
}

function map {

  # Check for the right number of args
  ExpectedNumArgs=5
  if [[ "$#" -ne "$ExpectedNumArgs" ]]; then
    echo "map function called with $# args; expected $ExpectedNumArgs."\
    "Quitting." >&2
    return 1
  fi

  # Assign the args
  ReadsToMap1=$1
  ReadsToMap2=$2
  LocalRef=$3
  OutFileStem=$4
  BamOnlyArg=$5

  # Some out files we'll produce
  InsertSizeCounts="$OutFileStem$InsertSizeCountsSuffix"
  Consensus="$OutFileStem"'_consensus_MinCov_'"$MinCov1"'_'"$MinCov2.fasta"
  BaseFreqs="$OutFileStem$BaseFreqsSuffix"
  BaseFreqsWHXB2="$OutFileStem$BaseFreqsWHXB2Suffix"
  ConsensusWcontigs="$OutFileStem"'_consensus_MinCov_'"$MinCov1"'_'"$MinCov2"'_wContigs.fasta'
  PreDedupBam="$OutFileStem$PreDeduplicationBamSuffix.bam"
  FinalOutBam="$OutFileStem.bam"
  DedupStats="$OutFileStem$DeduplicationStatsSuffix"

  # Check there's one seq in the ref file.
  NumSeqsInRefFile=$(grep -e '^>' "$LocalRef" | wc -l)
  if [ "$NumSeqsInRefFile" -eq 0 ]; then
    echo "Error: there are $NumSeqsInRefFile seqs in $LocalRef; there should"\
    "be exactly 1 for it to be used as a reference for mapping. Quitting." >&2
    return 1
  fi
  LocalRefName=$(awk '/^>/ {print substr($1,2)}' "$LocalRef")

  # Index the ref
  "$samtools" faidx "$LocalRef" && ls "$LocalRef".fai > /dev/null ||
  { echo 'Problem indexing the refererence with samtools. Quitting.' >&2 ; 
  return 1 ; }

  # Set the file name of the bam according to whether we will be deduplicating
  # it or not.
  if [[ "$deduplicate" == true ]]; then
    FinalConversionStepOut="$PreDedupBam"
  else
    FinalConversionStepOut="$FinalOutBam"
  fi


  # Map with the chosen mapper.
  if [[ "$mapper" == "smalt" ]]; then
    map_with_smalt "$ReadsToMap1" "$ReadsToMap2" "$LocalRef" \
    "$FinalConversionStepOut" ||
    { echo 'Problem mapping with smalt.' >&2 ; return 1 ; }
  elif [[ "$mapper" == "bowtie" ]]; then
    map_with_bowtie "$ReadsToMap1" "$ReadsToMap2" "$LocalRef" \
    "$FinalConversionStepOut" ||
    { echo 'Problem mapping with bowtie.' >&2 ; return 1 ; }
  elif [[ "$mapper" == "bwa" ]]; then
    map_with_bwa_mem "$ReadsToMap1" "$ReadsToMap2" "$LocalRef" \
    "$FinalConversionStepOut" ||
    { echo 'Problem mapping with bwa mem.' >&2 ; return 1 ; }
  else
    echo "Unrecognised value $mapper for the 'mapper' config file variable;"\
    "possible values are 'smalt', 'bowtie' or 'bwa'." >&2
    return 1
  fi

  # Deduplicate if desired
  if [[ "$deduplicate" == true ]]; then
    $DeduplicationCommand REMOVE_DUPLICATES=True I="$FinalConversionStepOut" \
    O="$FinalOutBam" M="$DedupStats" &&
    ls "$FinalOutBam" > /dev/null ||
    { echo "Problem running $DeduplicationCommand" >&2 ; return 1 ; }
  fi

  # Index the bam
  "$samtools" index "$FinalOutBam" ||
  { echo "Problem running $samtools index" >&2 ; return 1 ; }

  # Stop here if desired
  if [[ "$BamOnlyArg" == true ]]; then
    return 0
  fi

  # Check at least one read was mapped
  NumMappedReads=$(samtools view "$FinalOutBam" | wc -l)
  if [[ $NumMappedReads -eq 0 ]]; then
    echo "$FinalOutBam is empty - no reads were mapped!"
    return 1
  fi

  # Calculate the normalised insert size distribution.
  "$samtools" view "$FinalOutBam" | awk '{if ($9 > 0) print $9}' > "$InsertSizes1"
  InsertCount=$(wc -l "$InsertSizes1" | awk '{print $1}')
  if [[ $InsertCount -gt 0 ]]; then
    sort -n "$InsertSizes1" | uniq -c > "$InsertSizes2"
    awk '{print $2 "," $1 "," $1/'$InsertCount'}' "$InsertSizes2" > \
    "$InsertSizeCounts"
  else
    echo "Warning: no read in $FinalOutBam was identified as having"\
    "positive insert size. Unexpected. We'll skip making an insert size"\
    "distribution and continue."
  fi

  # Generate pileup
  echo 'Now calculating pileup - typically a slow step.'
  "$samtools" mpileup $mpileupOptions -f "$LocalRef" "$FinalOutBam" > \
  "$PileupFile" || { echo 'Failed to generate pileup.' >&2 ; return 1 ; }

  # Generate the base frequencies
  "$Code_AnalysePileup" "$PileupFile" "$LocalRef" > "$BaseFreqs" || \
  { echo 'Problem analysing the pileup.' >&2 ; return 1 ; }

  # Generate a version of the base freqs file with HXB2 coordinates, if desired.
  if [[ "$GiveHXB2coords" == "true" ]]; then
    HXB2file="$ThisDir/info/B.FR.83.HXB2_LAI_IIIB_BRU.K03455.fasta"
    if [[ ! -f "$HXB2file" ]]; then
      echo "Warning: the HXB2 sequence file, expected to be at $HXB2file, was"\
      "not found. We will not generate a version of the base frequency file"\
      "with HXB2 coordinates." >&2;
    else
      cat "$LocalRef" > "$RefWHXB2unaln"
      echo >> "$RefWHXB2unaln"
      cat "$HXB2file" >> "$RefWHXB2unaln"
      "$mafft" "$RefWHXB2unaln" > "$RefWHXB2aln" ||
      { echo "Problem running $mafft" >&2 ; return 1 ; }
      "$Code_MergeBaseFreqsAndCoords" "$BaseFreqs" --pairwise-aln \
      "$RefWHXB2aln" > "$BaseFreqsWHXB2" ||
      { echo "Problem running $Code_MergeBaseFreqsAndCoords" >&2 ; return 1 ; }
    fi
  fi

  # Call the consensuses
  "$Code_CallConsensus" "$BaseFreqs" "$MinCov1" "$MinCov2" "$MinBaseFrac" \
  --consensus-seq-name "$OutFileStem"'_consensus' --ref-seq-name "$LocalRefName" > \
  "$Consensus" || \
  { echo 'Problem calling the consensus.' >&2 ; return 1 ; }

  # Add the contigs to the alignment of the consensus and its reference.
  if [[ "$AlignContigsToConsensus" == "true" ]] && \
  [[ $NumHIVContigs -gt 0 ]]; then
    SwapContigsToTop=false
    AlignContigsToRefs "$Code_AlignToConsensus" '-S' "$RawContigFile2" \
    "$Consensus" "$ConsensusWcontigs" "$SwapContigsToTop" "$OldMafft" || \
    { echo 'Problem aligning the contigs to the consensus.' >&2 ; \
    return 1 ; }
  fi

}

# Designed for convenient use outside of shiver.
# Call it with three arguments: first the reads, second the reference sequence
# to which to map, thirdly a base name or stem for naming output files.
function MapUnpairedReadsStandAlone {

  # Options
  smalt='smalt'
  smaltIndexOptions="-k 15 -s 3"
  samtools='samtools'
  smaltMapOptions="-y 0.7"
  samtoolsReadFlags='-F 4'

  # temp files we'll create
  smaltIndex='temp_smaltRefIndex'
  MapOutConversion1='temp_MapOutStep1'
  MapOutConversion2='temp_MapOutStep2'
  MapOutAsSam='temp_MapOut.sam'

  # Check for the right number of args
  ExpectedNumArgs=3
  if [[ "$#" -ne "$ExpectedNumArgs" ]]; then
    echo "MapUnpairedReadsStandAlone function called with $# args; expected"\
    "$ExpectedNumArgs. Quitting." >&2
    return 1
  fi

  # Assign the args
  ReadsToMap=$1
  LocalRef=$2
  OutFileStem=$3

  # Check there's one seq in the ref file.
  NumSeqsInRefFile=$(grep -e '^>' "$LocalRef" | wc -l)
  if [ "$NumSeqsInRefFile" -eq 0 ]; then
    echo "Error: there are $NumSeqsInRefFile seqs in $LocalRef; there should"\
    "be exactly 1 for it to be used as a reference for mapping. Quitting." >&2
    return 1
  fi
  LocalRefName=$(awk '/^>/ {print substr($1,2)}' "$LocalRef")

  # Index the ref
  "$smalt" index $smaltIndexOptions "$smaltIndex" "$LocalRef" ||
  { echo 'Problem indexing the refererence with smalt. Quitting.' >&2 ;
  return 1 ; }
  "$samtools" faidx "$LocalRef" ||
  { echo 'Problem indexing the refererence with samtools. Quitting.' >&2 ; 
  return 1 ; }

  # Do the mapping!
  echo 'Now mapping - typically a slow step.'
  "$smalt" map $smaltMapOptions -o "$MapOutAsSam" "$smaltIndex" \
  "$ReadsToMap" || { echo 'Smalt mapping failed.' >&2 ; return 1 ; }

  # Convert that sam file into a bam file. Thanks Nick Croucher!
  "$samtools" view -bS $samtoolsReadFlags -t "$LocalRef".fai -o \
  "$MapOutConversion1".bam "$MapOutAsSam" &&
  "$samtools" sort "$MapOutConversion1".bam -o "$OutFileStem".bam -T "$SamtoolsSortFile" &&
  "$samtools" index "$FinalOutBam" || \
  { echo 'Failed to convert from sam to bam format.' >&2 ; return 1 ; }

}

function GetHIVcontigs {

  # Check for the right number of args and assign them.
  ExpectedNumArgs=4
  if [[ "$#" -ne "$ExpectedNumArgs" ]]; then
    echo "GetHIVcontigs function called with $# args; expected"\
    "$ExpectedNumArgs. Quitting." >&2
    return 1
  fi
  ContigFile=$1
  LongContigs=$2
  BlastFile=$3
  HIVcontigsFile=$4

  # Check that there are some contigs, that their IDs are unique, and that their
  # IDs don't contain commas.
  ContigNames=$(awk '/^>/ {print substr($1,2)}' "$ContigFile" | sort)
  NumContigs=$(echo "$ContigNames" | wc -w)
  if [[ $NumContigs -eq 0 ]]; then
    echo "$ContigFile contains no sequences." >&2
    return 1;
  fi
  NumUniqueIDs=$(printf '%s\n' $ContigNames | uniq | wc -l)
  if [[ $NumUniqueIDs -ne $NumContigs ]]; then
    echo "$ContigFile contains some identically named sequences. Rename"\
    'these and try again.' >&2
    return 1;
  fi
  if [[ "$ContigNames" == *","* ]]; then
    echo "Contig names must not contain commas." >&2
    return 1
  fi

  # Create a new file of contigs keeping only those that are long enough.
  "$Code_FindSeqsInFasta" "$ContigFile" "" --match-start --min-length \
  "$MinContigLength" > "$LongContigs" || { echo "Problem removing short"\
  "contigs from $ContigFile." >&2 ; return 1 ; }

  # Blast the contigs. Keep only hits for which (length * fractional identity)
  # is longer than the minimum contig length.
  "$BlastNcommand" -query "$LongContigs" -db "$BlastDatabase" -max_target_seqs \
  1 -outfmt '10 qseqid sseqid evalue pident qlen qstart qend sstart send' \
  -word_size "$BlastWordSize" | awk -F, \
  '($4/100 * ($7-$6+1))>='"$MinContigLength" > "$BlastFile" || 
  { echo "Problem blasting $LongContigs." >&2 ; return 1 ; }

  # If there are no blast hits, nothing needs doing. Exit.
  NumBlastHits=$(wc -l "$BlastFile" | awk '{print $1}')
  if [ "$NumBlastHits" -eq 0 ]; then
    echo "No contig in $LongContigs has a blast hit against any of the"\
    "references used to create the shiver initialisation directory (i.e. this"\
    "sample is presumably pure contamination)."
    return 1
  fi

  # Extract those contigs that have a blast hit.
  HIVcontigNames=$(awk -F, '{print $1}' "$BlastFile" | sort | uniq)
  "$Code_FindSeqsInFasta" "$LongContigs" $HIVcontigNames > \
  "$HIVcontigsFile" || { echo 'Problem extracting the HIV contigs using the'\
  'blast results.' >&2 ; return 1 ; }

}


function CheckConfig {

  # Check for the right number of args and assign them.
  ExpectedNumArgs=4
  if [[ "$#" -ne "$ExpectedNumArgs" ]]; then
    echo "CheckConfig function called with $# args; expected"\
    "$ExpectedNumArgs. Quitting." >&2
    return 1
  fi
  ConfigFile=$1
  CheckForInit=$2
  CheckForAligningContigs=$3
  CheckForMapping=$4
  source "$ConfigFile"

  # Coerce boolean args into true bools if needed
  if [[ "$CheckForInit" == "true" ]]; then
    CheckForInit=true
  else
    CheckForInit=false
  fi
  if [[ "$CheckForAligningContigs" == "true" ]]; then
    CheckForAligningContigs=true
  else
    CheckForAligningContigs=false
  fi
  if [[ "$CheckForMapping" == "true" ]]; then
    CheckForMapping=true
  else
    CheckForMapping=false
  fi

  # Check 0 < MaxContigGappiness < 1, and 0 < MinContigHitFrac < 1
  NonNegativeRegex='^[0-9]+([.][0-9]+)?$'
  if [[ "$MaxContigGappiness" =~ $NonNegativeRegex ]] && \
  (( $(echo "$MaxContigGappiness > 0" | bc -l) )) && \
  (( $(echo "$MaxContigGappiness < 1" | bc -l) )); then
    :
  else
    echo "The 'MaxContigGappiness' variable in the config file should be" \
    "greater than 0 and less than 1." >&2
    return 1
  fi
  if [[ "$MinContigHitFrac" =~ $NonNegativeRegex ]] && \
  (( $(echo "$MinContigHitFrac > 0" | bc -l) )) && \
  (( $(echo "$MinContigHitFrac < 1" | bc -l) )); then
    :
  else
    echo "The 'MinContigHitFrac' variable in the contig file should be a" \
    "number greater than 0 and less than 1." >&2
    return 1
  fi

  # Check BlastDBcommand works
  if $CheckForInit || $CheckForMapping; then
    "$BlastDBcommand" -help &> /dev/null || { echo "Error running" \
    "'$BlastDBcommand -help'. Are you sure that blast is installed, and that you"\
    "chose the right value for the config file variable 'BlastDBcommand'?" >&2; \
    return 1; }
  fi

  # Check BlastNcommand works
  if $CheckForAligningContigs || $CheckForMapping; then
    "$BlastNcommand" -help &> /dev/null || { echo "Error running" \
    "'$BlastNcommand -help'. Are you sure that blast is installed, and that"\
    "you chose the right value for the config file variable 'BlastNcommand'?"\
    >&2; return 1; }
  fi

  # Check samtools works. NB older versions don't have the help option, so try
  # view the test sam file.
  if $CheckForMapping; then
    "$samtools" help &> /dev/null ||
    "$samtools" view "$ToolsDir"/test.sam -S &> /dev/null || { echo "Could" \
    'not run either'
    echo "$samtools" help
    echo 'or'
    echo "$samtools" view "$ToolsDir"/test.sam -S
    echo 'Are you sure that samtools is installed, and that you chose the' \
    "right value for the config file variable 'samtools'? Quitting." >&2; \
    return 1; }
  fi

  # Check mafft works
  if $CheckForAligningContigs || $CheckForMapping; then
    echo -e ">seq1\naa\n>seq2\naa" | "$mafft" - &> /dev/null || { echo "Error" \
    "running '$mafft'. Are you sure that mafft is installed, and that you"\
    "chose the right value for the config file variable 'mafft'?" >&2; \
    return 1; }
  fi

  # Check boolean variables are either true or false.
  if [[ "$TrimReadsForAdaptersAndQual" != "true" ]] && \
  [[ "$TrimReadsForAdaptersAndQual" != "false" ]]; then
    echo "The 'TrimReadsForAdaptersAndQual' variable in the config file should"\
    "be either true or false."
    return 1
  fi
  if [[ "$TrimReadsForPrimers" != "true" ]] && \
  [[ "$TrimReadsForPrimers" != "false" ]]; then
    echo "The 'TrimReadsForPrimers' variable in the config file should"\
    "be either true or false."
    return 1
  fi
  if [[ "$CleanReads" != "true" ]] && [[ "$CleanReads" != "false" ]]; then
    echo "The 'CleanReads' variable in the config file should"\
    "be either true or false."
    return 1
  fi
  if [[ "$remap" != "true" ]] && [[ "$remap" != "false" ]]; then
    echo "The 'remap' variable in the config file should"\
    "be either true or false."
    return 1
  fi
  if [[ "$MapContaminantReads" != "true" ]] && \
  [[ "$MapContaminantReads" != "false" ]]; then
    echo "The 'MapContaminantReads' variable in the config file should"\
    "be either true or false."
    return 1
  fi

  # Some checks only needed if we're mapping:
  if $CheckForMapping; then

    # Check fastaq works, if needed
    if [[ "$TrimReadsForPrimers" == "true" ]]; then
      "$fastaq" version &> /dev/null || { echo "Error running" \
      "'$fastaq version'. Are you sure that fastaq is installed, and that you"\
      "chose the right value for the config file variable 'fastaq'?" >&2; \
      return 1; }
    fi

    # Check trimmomatic works, if needed
    if [[ "$TrimReadsForAdaptersAndQual" == "true" ]]; then
      $trimmomatic -version &> /dev/null || { echo "Error running" \
      "'$trimmomatic -version'. Are you sure that trimmomatic is installed,"\
      "and that you chose the right value for the config file variable" \
      "'trimmomatic'?" >&2; return 1; }
    fi

    # Test the mapper is one we support and can run.
    if [[ "$mapper" == "smalt" ]]; then
      "$smalt" version &> /dev/null || { echo "Error running" \
      "'$smalt version'. Are you sure that smalt is installed, and that you"\
      "chose the right value for the config file variable 'mapper'?" >&2; \
      return 1; }
    elif [[ "$mapper" == "bowtie" ]]; then
      "$bowtie2" --help &> /dev/null || { echo "Error running" \
      "$bowtie2 --help. Are you sure that bowtie2 is installed, and"\
      "that you chose the right value for the config file variable" \
      "'mapper'?" >&2; return 1; }
    elif [[ "$mapper" == "bwa" ]]; then
      # bwa doesn't seem to have any kind of 'help' or 'version' command that we
      # can call to test it works.
      :
    else
      echo "Unrecognised value $mapper for the 'mapper' config file variable;"\
      "possible values are 'smalt', 'bowtie' or 'bwa'." >&2
      return 1
    fi
  fi

  # Check positive ints are positive ints
  NonNegativeIntRegex='^[0-9]+$'
  if ! [[ "$NumThreadsTrimmomatic" =~ $NonNegativeIntRegex ]] || \
  [[ "$NumThreadsTrimmomatic" -lt 1 ]]; then
    echo "The 'NumThreadsTrimmomatic' variable in the config file should be an"\
    "integer greater than 0." >&2
    return 1
  fi
  if ! [[ "$MinCov1" =~ $NonNegativeIntRegex ]] || \
  [[ "$MinCov1" -lt 1 ]]; then
    echo "The 'MinCov1' variable in the config file should be an"\
    "integer greater than 0." >&2
    return 1
  fi
  if ! [[ "$MinCov2" =~ $NonNegativeIntRegex ]] || \
  [[ "$MinCov2" -lt 1 ]]; then
    echo "The 'MinCov2' variable in the config file should be an"\
    "integer greater than 0." >&2
    return 1
  fi
  if ! [[ "$MinContigLength" =~ $NonNegativeIntRegex ]] || \
  [[ "$MinContigLength" -lt 1 ]]; then
    echo "The 'MinContigLength' variable in the config file should be an"\
    "integer greater than 0." >&2
    return 1
  fi

  # Check MinCov2 >= MinCov1
  if [[ "$MinCov2" -lt "$MinCov1" ]]; then
    echo "The 'MinCov2' variable in the config file should be greater than or"\
    "equal to MinCov1." >&2
    return 1
  fi

  # Check MinBaseFrac <= 1
  FloatRegex='^-?[0-9]+([.][0-9]+)?$'
  if [[ "$MinBaseFrac" =~ $FloatRegex ]] && \
  (( $(echo "$MinBaseFrac <= 1" | bc -l) )); then
    :
  else
    echo "The 'MinBaseFrac' variable in the config file should be equal to or" \
    "less than 1." >&2
    return 1
  fi

}
