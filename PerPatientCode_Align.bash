#!/usr/bin/env bash

# Arguments for this script:
# 1) a fasta file of contigs (output from processing the short reads with an
# assembly program),
# 2) an alignment of references,
# 3) a basename used for naming the output from this script (a sensible choice
# might be the contig file name minus its path and extension).
# If this script completes successfully it will produce a .blast file (detailing
# the contigs that have blast hits). If the .blast file is not empty it will
# produce a fasta file of those contigs with hits and another fasta file of
# these contigs aligned to the refereces; if cutting and/or reversing of contigs
# is necessary, two more fasta files are produced - the cut/reversed contigs on
# their own and also aligned to references (i.e. there will be two files of the
# contigs on their own and two files of the contigs aligned to references).

set -u

################################################################################
# INITIALISATION

# To be created by a 'pipeline prep' one-off step?
BlastDataBase="$HOME/JobInputs/blastDB_HIV1_COM_2012_genome_DNA.db"

# Some code we'll need:
FindSeqsInFasta="$HOME/AnalysisCode/FindSeqsInFasta.py"
ContigCutter="$HOME/AnalysisCode/sort.blast.output.standalone.R"
FastaFileComparer="$HOME/AnalysisCode/CheckFastaFileEquality.py"

# Check for the right number of arguments. Assign them to variables.
NumArgsExpected=3
if [ "$#" -ne "$NumArgsExpected" ]; then
  echo "$#" 'arguments specified;' "$NumArgsExpected" 'expected. Quitting' >&2
  exit 1
fi
ContigFile="$1"
RefAlignment="$2"
OutFileBasename="$3"

# Check that certain input files exist
function CheckFilesExist {
  for argument in "$@"; do
    if [ ! -f "$argument" ]; then
      echo "$argument" 'does not exist. Quitting.' >&2
      exit 1
    fi
  done
}
CheckFilesExist "$ContigFile" "$RefAlignment" 

# The names for output files we'll produce.
BlastFile="$OutFileBasename".blast
RawContigFile="$OutFileBasename"_hiv.fasta
CutContigFile="$OutFileBasename"_hiv_cut.fasta
RawContigAlignment="$OutFileBasename"_raw_wRefs.fasta
CutContigAlignment="$OutFileBasename"_cut_wRefs.fasta
################################################################################


################################################################################
# THE MAIN PROGRAM

# Blast the contigs
blastn -query "$ContigFile" -db "$BlastDataBase" -out "$BlastFile" \
-max_target_seqs 1 -outfmt \
'10 qacc sacc sseqid evalue pident qstart qend sstart send' || \
{ echo 'Problem blasting' "$ContigFile"'. Quitting.' >&2 ; exit 1 ; }

# If there are no blast hits, nothing needs doing. Exit.
NumBlastHits=$(wc -l "$BlastFile" | awk '{print $1}')
if [ "$NumBlastHits" -eq 0 ]; then
  exit 0
fi

# Extract those contigs that have a blast hit...
awk -F, '{print $1}' "$BlastFile" | sort | uniq | xargs \
"$FindSeqsInFasta" "$ContigFile" > "$RawContigFile" || \
{ echo 'Problem extracting the HIV contigs. Quitting.' >&2 ; exit 1 ; }

# ...and align them to the refs.
mafft --quiet --add "$RawContigFile" "$RefAlignment" \
> "$RawContigAlignment" || \
{ echo 'Problem aligning' "$ContigFile"'. Quitting.' >&2 ; exit 1 ; }

# Run the contig cutting & flipping code.
if Rscript "$ContigCutter" "$BlastFile" "$ContigFile" "$CutContigFile"; then
  # The contig cutting & flipping code worked.

  "$FastaFileComparer" "$RawContigFile" "$CutContigFile"
  ComparisonExitStatus=$?

  if [ $ComparisonExitStatus -eq 0 ]; then
    # The contig cutting & flipping code left the contigs untouched. Do nothing.
    echo 'The contig cutting & flipping code found nothing to do.'

  elif [ $ComparisonExitStatus -eq 111 ]; then
    # The contig cutting & flipping code modified the contigs. Align them.
    mafft --quiet --add "$CutContigFile" "$RefAlignment" \
    > "$CutContigAlignment" || \
    { echo 'Problem aligning' "$CutContigFile"'. Quitting.' >&2 ; exit 1 ; }

  else
    echo 'Problem running' "$FastaFileComparer"'. Quitting.' >&2
    exit 1
  fi  

else
  echo 'Error encountered running' "$ContigCutter"
fi
################################################################################
