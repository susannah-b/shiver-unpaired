#!/usr/bin/env bash

ThisDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ToolsDir="$ThisDir"/tools
Code_AnalysePileup="$ToolsDir/AnalysePileup.py"
Code_CallConsensus="$ToolsDir/CallConsensus.py"
Code_CheckFastaFileEquality="$ToolsDir/CheckFastaFileEquality.py"
Code_ConstructRef="$ToolsDir/ConstructBestRef.py"
Code_CorrectContigs="$ToolsDir/CorrectContigs.py"
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

function AlignContigs {

  ContigFile=$1
  OutputAlignedContigFile=$2
  OldMafftArg=$3
  ContigNames=$(awk '/^>/ {print substr($1,2)}' "$ContigFile")

  mafft --quiet --add "$ContigFile" "$RefAlignment" > "$TempContigAlignment1" \
  || { echo 'Problem aligning' "$ContigFile"'. Quitting.' >&2 ; exit 1 ; }
  MaxContigGappiness1=$("$Code_ConstructRef" -S1 "$TempContigAlignment1" \
  $ContigNames | sort -nrk2,2 | head -1 | awk '{print $2}')
  BestContigAlignment="$TempContigAlignment1"

  # If the function was called with OldMafftArg=true, don't bother trying the
  # addfragments option. Otherwise, try. If it doesn't work, set OldMafft=true
  # ready for the next call of this function; if it does work, use the least
  # gappy alignment.
  if ! $OldMafftArg; then
    mafft --quiet --addfragments "$ContigFile" "$RefAlignment" \
    > "$TempContigAlignment2" || \
    { echo "Warning: it looks like you're running an old version of mafft: the" \
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
      fi
    fi
  fi

  # Swap the contigs from after the references to before them, for easier visual
  # inspection.
  NumRefs=$(grep -e '^>' "$RefAlignment" | wc -l)
  ContigsStartLine=$(awk '/^>/ {N++; if (N=='$((NumRefs+1))') {print NR; exit}}' \
  "$BestContigAlignment")
  tail -n +"$ContigsStartLine" "$BestContigAlignment" > "$OutputAlignedContigFile"
  head -n "$((ContigsStartLine-1))" "$BestContigAlignment" >> \
  "$OutputAlignedContigFile"

}
