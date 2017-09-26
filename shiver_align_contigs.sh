#!/usr/bin/env bash

set -u
set -o pipefail

UsageInstructions=$(echo '
Arguments for this script:
(1) the initialisation directory you created using the shiver_init.sh command;
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
  echo "$InitDir does not exist or is not a directory. Quitting." >&2
  exit 1
fi
InitDir=$(cd "$InitDir"; pwd)

BlastDatabase="$InitDir/ExistingRefsBlastDatabase"
RefAlignment="$InitDir/ExistingRefAlignment.fasta"

# Source the shiver funcs, check files exist, source the config file, check it.
ThisDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source "$ThisDir"/'shiver_funcs.sh'
CheckFilesExist "$ContigFile" "$RefAlignment"
CheckConfig "$ConfigFile" || \
{ echo "Problem with $ConfigFile. Quitting." >&2 ; exit 1 ; }

# Out files we'll make
LongContigs="$SID$LongEnoughContigsSuffix"
BlastFile="$SID$BlastSuffix"
RawContigAlignment="$SID"'_raw_wRefs.fasta'
CutContigAlignment="$SID"'_cut_wRefs.fasta'

# Extract just the HIV contigs (those that blast to the refs) and put them in
# $RawContigFile1
GetHIVcontigs "$ContigFile" "$LongContigs" "$BlastFile" "$RawContigFile1" || \
{ echo "Problem encountered while checking the contigs in $ContigFile and"\
" extracting those thought to be HIV. Quitting." >&2 ; exit 1 ; }

# Align the HIV contigs.
OldMafft=false
SwapContigsToTop=true
AlignContigsToRefs "$mafft" '--quiet' "$RawContigFile1" "$RefAlignment" \
"$RawContigAlignment" "$SwapContigsToTop" "$OldMafft" || \
{ echo 'Problem aligning the raw contigs to refs. Quitting.' >&2 ; exit 1 ; }

# Run the contig cutting & flipping code
"$Code_CorrectContigs" "$BlastFile" -C "$ContigFile" -O "$CutContigFile" || \
{ echo "Problem encountered running $Code_CorrectContigs. Quitting." >&2 ; \
exit 1; }

# If the contigs needed cutting and/or flipping: align the modified contigs.
if [[ -f "$CutContigFile" ]]; then
  AlignContigsToRefs "$mafft" '--quiet' "$CutContigFile" "$RefAlignment" \
  "$CutContigAlignment" "$SwapContigsToTop" "$OldMafft" || \
  { echo 'Problem aligning the cut/modified contigs to refs. Quitting.' >&2 ;
  exit 1 ; }
fi

