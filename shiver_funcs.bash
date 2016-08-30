#!/usr/bin/env bash

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
  { echo 'Problem aligning' "$ContigFile"'. Quitting.' >&2 ; exit 1 ; }
  MaxContigGappiness1=$("$Code_ConstructRef" -S1 "$TempContigAlignment1" \
  $ContigNames | sort -nrk2,2 | head -1 | awk '{print $2}')
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
      $ContigNames | sort -nrk2,2 | head -1 | awk '{print $2}')
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
  # Forward reads (1) or reverse reads (2)
  OneOrTwo=$2

  # Check all read seq ids end in /1 or /2 as needed.
  suffix=$(awk '{if ((NR-1)%4==0) print substr($1,length($1)-1,
  length($1))}' "$ReadFile" | sort | uniq)
  if [[ "$suffix" != '/'"$OneOrTwo" ]]; then
    echo "Found at least one read in $ReadFile whose sequence ID does not end"\
    'in "\'"$OneOrTwo"'". Quitting.'
    exit 1
  fi

  # Check none of the lines with read IDs contain tabs
  NumNameLinesWithTabs=$(awk '{if ((NR-1)%4==0 && gsub("\t","\t",$0) > 0)
  print}' "$ReadFile" | wc -l)
  if [[ $NumNameLinesWithTabs -ne 0 ]]; then
    echo "The following lines in $ReadFile contain tabs:"
    awk '{if ((NR-1)%4==0 && gsub("\t","\t",$0) > 0) print}' "$ReadFile"
    echo 'To remove contaminant reads, we require there to be no tabs in the' \
    'sequence ID lines of fastq files. Quitting.'
    exit 1
  fi

}
