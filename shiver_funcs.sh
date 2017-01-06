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

# For quitting if files don't exist.
function CheckFilesExist {
  for argument in "$@"; do
    if [ ! -f "$argument" ]; then
      echo "$argument" 'does not exist. Quitting.' >&2
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
    if [[ $NumLinesInAln -gt 0 ]]; then 
      MaxContigGappiness2=$("$Code_ConstructRef" -S1 "$TempContigAlignment2" \
      $ContigNames | sort -nrk2,2 | head -1 | awk '{print $2}') || { echo \
      "Problem analysing $TempContigAlignment1 (i.e. the output from aligning "\
      "$ContigFile and $ThisRefAlignment) with $Code_ConstructRef." \
      >&2 ; return 1 ; }
      if (( $(echo "$MaxContigGappiness2 < $MaxContigGappiness1" | bc -l) )); 
      then
        BestContigAlignment="$TempContigAlignment2"
        echo 'Info: the --addfragments mafft option performed better than'\
        '--add. Using the former.'
      else
        echo 'Info: the --addfragments mafft option performed no better than'\
        '--add. Using the latter.'
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

function map {

  LocalRef=$1
  LocalRefName=$2
  OutFileStem=$3

  InsertSizeCounts="$OutFileStem$InsertSizeCountsSuffix"
  Consensus="$OutFileStem"'_consensus_MinCov_'"$MinCov1"'_'"$MinCov2.fasta"
  BaseFreqs="$OutFileStem$BaseFreqsSuffix"
  ConsensusWcontigs="$OutFileStem"'_consensus_MinCov_'"$MinCov1"'_'"$MinCov2"'_wContigs.fasta'

  # Do the mapping!
  echo 'Now mapping - typically a slow step.'
  "$smalt" map $smaltMapOptions -o "$MapOutAsSam" "$smaltIndex" \
  "$cleaned1reads" "$cleaned2reads" || \
  { echo 'Smalt mapping failed.' >&2 ; return 1 ; }

  # Convert that sam file into a bam file. Thanks Nick Croucher!
  "$samtools" view -bS $samtoolsReadFlags -t "$LocalRef".fai -o \
  "$MapOutConversion1".bam "$MapOutAsSam" &&
  "$samtools" sort -n "$MapOutConversion1".bam -o "$MapOutConversion2".bam -T "$SamtoolsSortFile" &&
  "$samtools" fixmate "$MapOutConversion2".bam "$MapOutConversion3".bam &&
  "$samtools" sort "$MapOutConversion3".bam -o "$OutFileStem".bam -T "$SamtoolsSortFile" &&
  "$samtools" index "$OutFileStem.bam" || \
  { echo 'Failed to convert from sam to bam format.' >&2 ; return 1 ; }

  # Check at least one read was mapped
  NumMappedReads=$(samtools view "$OutFileStem.bam" | wc -l)
  if [[ $NumMappedReads -eq 0 ]]; then
    echo "$OutFileStem.bam is empty - no reads were mapped!"
    return 1
  fi

  # Calculate the normalised insert size distribution.
  "$samtools" view "$OutFileStem.bam" | awk '{if ($9 > 0) print $9}' > "$InsertSizes1"
  InsertCount=$(wc -l "$InsertSizes1" | awk '{print $1}')
  if [[ $InsertCount -gt 0 ]]; then
    sort -n "$InsertSizes1" | uniq -c > "$InsertSizes2"
    awk '{print $2 "," $1 "," $1/'$InsertCount'}' "$InsertSizes2" > \
    "$InsertSizeCounts"
  else
    echo "Warning: no read in $OutFileStem.bam was identified as having"\
    "positive insert size. Unexpected. We'll skip making an insert size"\
    "distribution and continue."
  fi

  # Generate pileup
  echo 'Now calculating pileup - typically a slow step.'
  "$samtools" mpileup $mpileupOptions -f "$LocalRef" "$OutFileStem.bam" > \
  "$PileupFile" || { echo 'Failed to generate pileup.' >&2 ; return 1 ; }

  # Generate base frequencies and the consensuses
  "$Code_AnalysePileup" "$PileupFile" "$LocalRef" > "$BaseFreqs" || \
  { echo 'Problem analysing the pileup.' >&2 ; return 1 ; }
  "$Code_CallConsensus" "$BaseFreqs" "$MinCov1" "$MinCov2" "$MinBaseFrac" \
  --consensus-seq-name "$OutFileStem"'_consensus' --ref-seq-name "$LocalRefName" > \
  "$Consensus" || \
  { echo 'Problem calling the consensus.' >&2 ; return 1 ; }

  if [[ $NumHIVContigs -gt 0 ]]; then
    # Add the contigs to the alignment of the consensus and its reference.
    SwapContigsToTop=false
    AlignContigsToRefs "$Code_AlignToConsensus" '-S' "$RawContigFile2" \
    "$Consensus" "$ConsensusWcontigs" "$SwapContigsToTop" "$OldMafft" || \
    { echo 'Problem aligning the contigs to the consensus.' >&2 ; \
    return 1 ; }
  fi

}


function CheckConfig {

  ConfigFile="$1"
  source "$ConfigFile"

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
  "$BlastDBcommand" -help &> /dev/null || { echo "Error running" \
  "'$BlastDBcommand -help'. Are you sure that blast is installed, and that you"\
  "chose the right value for the config file variable 'BlastDBcommand'?" >&2; \
  return 1; }

  # Check BlastNcommand works
  "$BlastNcommand" -help &> /dev/null || { echo "Error running" \
  "'$BlastNcommand -help'. Are you sure that blast is installed, and that you"\
  "chose the right value for the config file variable 'BlastNcommand'?" >&2; \
  return 1; }

  # Check smalt works
  "$smalt" version &> /dev/null || { echo "Error running" \
  "'$smalt version'. Are you sure that smalt is installed, and that you"\
  "chose the right value for the config file variable 'smalt'?" >&2; \
  return 1; }

  # Check samtools works
  "$samtools" help &> /dev/null || { echo "Error running" \
  "'$samtools help'. Are you sure that samtools is installed, and that you"\
  "chose the right value for the config file variable 'samtools'?" >&2; \
  return 1; }

  # Check mafft works
  echo -e ">seq1\naa\n>seq2\naa" | "$mafft" - &> /dev/null || { echo "Error" \
  "running '$mafft'. Are you sure that mafft is installed, and that you"\
  "chose the right value for the config file variable 'mafft'?" >&2; \
  return 1; }

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

  # Check fastaq works, if needed
  if [[ "$TrimReadsForPrimers" == "true" ]]; then
    "$fastaq" version &> /dev/null || { echo "Error running" \
    "'$fastaq version'. Are you sure that fastaq is installed, and that you"\
    "chose the right value for the config file variable 'fastaq'?" >&2; \
    return 1; }
  fi

  # Check trimmomatic works, if needed
  if [[ "$TrimReadsForAdaptersAndQual" == "true" ]]; then
    "$trimmomatic" -version &> /dev/null || { echo "Error running" \
    "'$trimmomatic -version'. Are you sure that trimmomatic is installed, and"\
    "that you chose the right value for the config file variable" \
    "'trimmomatic'?" >&2; return 1; }
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
