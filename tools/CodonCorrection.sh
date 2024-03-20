#!/usr/bin/env bash

# CodonCorrectionInit.sh must be run once before using this script

set -u
set -o pipefail
set -e

UsageInstructions="Arguments for this script:
(1) A sample file containing the sequence to be corrected as the first sequence;
(2) A file of gene coordinates - NB this requires specific formatting/order (see note in code)
(3) The filepath for the virulign binaries
(4) The export options for virulign separated by ',' and enclosed by double quotes. One of more of Mutations, MutationTable, \
Nucleotides, or Amino Acids, or use 'all' to generate all four.
(5) The chosen genes to analyse for the sample separated by ',' and enclosed by double quotes. \
One or more of gag, pol, env, vpr, vif, vpu, or nef, or use 'all' to analyse all.
(6) The init directory.
(7) The directory where output files will be created.
"
# Write help text to explain the different export options

# Source config - improve to not rely on specific filepath
source "$HOME/shiver/config.sh"

### Check arguments
ExpectedArgs=7
if [ $# -ne "$ExpectedArgs" ]; then
  echo "$UsageInstructions"
  echo "Incorrect number of arguments provided, expected:" "$ExpectedArgs". "Quitting." >&2
  exit 1
else
  SampleFile=$1
  GeneCoordInfo=$2
  VirulignLocation=$3
  VirulignOptions=$4
  GenesToAnalyse=$5
  InitDir=$6
  OutputDir=$7
fi

# Colours
RED="\033[91m"
BLUE="\033[94m"
YELLOW="\033[93m"
GREEN="\033[92m"
GREY="\033[37m"
END="\033[0m"


# Other variables
ReferenceFile="$InitDir/CC_References.fasta"

# Check init script has been run
BLASTnDB_Path="$InitDir/BLASTnDB"
if [[ ! -d "$BLASTnDB_Path" ]]; then
  echo "No BLASTn database found. Please run CodonCorrectionInit.sh before running this script. Quitting." >&2
  exit 1
fi

# Check if file exists
# Call as a function and do for all
if ! [ -f "$SampleFile" ]; then
  echo "The provided sample file '$SampleFile' does not exist. Check the file path is correct."
  exit 1
# also check other files
fi

# Check sample and reference are .fasta files
if [[ "$SampleFile" != *.fasta ]]; then 
  echo "Sample file $SampleFile is not a fasta file. Quitting." >&2
  exit 1
fi
# Check reference file is a .fasta
if [[ "$ReferenceFile" != *.fasta ]]; then 
  echo "Reference file $ReferenceFile is not a fasta file. Quitting." >&2
  exit 1
fi

### Check the sample and reference files
# Count sequence number in reference and sample
Ref_SeqNumber=$(grep '^>' "$ReferenceFile" | wc -l | awk '{$1=$1};1')
if [[ "$Ref_SeqNumber" == 0 ]]; then
  echo "Reference file $ReferenceFile contains no sequences. Quitting." >&2
  exit 1
fi
Sample_SeqNumber=$(grep '^>' "$SampleFile" | wc -l | awk '{$1=$1};1')
if [[ "$Sample_SeqNumber" == 0 ]]; then
  echo "Reference file $ReferenceFile contains no sequences. Quitting." >&2
  exit 1
fi
# Determine whether to use mafft pairwise or AlignMoreSeqs - perhaps implement more rigourous checks later
if [[ "$Sample_SeqNumber" == 1 ]]; then
  SingleSequence=true # Will be funneled into mafft pairwise alignment
elif [[ "$Sample_SeqNumber" == 2 ]]; then
  SingleSequence=false # In this case expects remap consensus files (=2 sequences) and will do AlignMoreSeqs. Can develop this later for other file formats if needed
else 
  echo "Unexpected number of sequences in the sample file. Expected 1 or 2. Quitting." >&2
  exit 1
fi

# Check GeneCoordInfo
if [[ "$GeneCoordInfo" != *.txt ]]; then 
  echo "Gene coordinates file in unexpected format. Please provide as a .txt file. Quitting." >&2
  exit 1
fi
# Check header is as expected
# Not strictly necessary for function but ensures the gene info is formatted in the correct way. But maybe overly strict.
GC_FirstLine=$(head -n 1 "$GeneCoordInfo")
GeneCoord_Header="Sequence_name, gag_start, gag_end, pol_start, pol_end, vif_start, vif_end, vpr_start, \
vpr_end, vpu_start, vpu_end, env_start, env_end, nef_start, nef_end"
if [[ "$GC_FirstLine" != *"$GeneCoord_Header"* ]]; then
  echo -e "Expected the gene coordinates file to have the header: '$GeneCoord_Header'\nPlease supply the gene coordinate \
  data in the correct format. Quitting" >&2
  exit 1
fi

# Check VIRULIGN exists
if [[ ! -f "$VirulignLocation" ]]; then
  echo "The specified VIRULIGN location '$VirulignLocation' does not exist. Check file path is correct. Quitting." >&2
  exit 1
fi

# Check VIRULIGN
"$VirulignLocation" &> /dev/null || echo "Error running VIRULIGN. Check the program is installed correctly." >&2

# cd to $OutputDir
if [[ -d "$OutputDir" ]]; then
  echo "$OutputDir exists already; quitting to prevent overwriting." >&2
  exit 1
fi
mkdir -p "$OutputDir" && cd "$OutputDir" ||
{ echo "Could not mkdir then cd to $OutputDir. Quitting." >&2 ; exit 1 ; }

# IMPORTANT: Gene Coords must be set up in the correct format. My scripts use gene fragments to determine gene coordinates
# from complete genomes and make a specific .txt, however I haven't integrated this into the CodonCorrection script. 
# If just providing a set of references with this tool, then perhaps parts could be simplified to reflect this

### Other variables
# Set up for $shiver
AlignMoreSeqsTool=~/shiver/tools/AlignMoreSeqsToPairWithMissingCoverage.py
# Extracting the sample ID from the fasta
SequenceName_shell=$(head -n 1 "$SampleFile" | sed 's/>//')

# Check config and other shiver scripts for variables to work from if needed

# Add red to all the error messages, yellow to warnings, blue for standard messages, green for results

###################################################################################

### Determine which genes to correct

GAG=false
POL=false
VIF=false
VPR=false
VPU=false
ENV=false
NEF=false

genes=("GAG" "POL" "VIF" "VPR" "VPU" "ENV" "NEF")

# Switch on all genes if command set to "all"
if [[ "$GenesToAnalyse" == "all" ]]; then
  GAG=true
  POL=true
  VIF=true
  VPR=true
  VPU=true
  ENV=true
  NEF=true
else
  # Switch on the corresponding gene variable for each listed gene
  IFS=',' read -ra listed_genes <<< "$GenesToAnalyse"
  for gene_name in "${listed_genes[@]}"; do
      # Convert gene name to uppercase
      gene_upper=$(echo "$gene_name" | tr '[:lower:]' '[:upper:]')

      # Check if the gene name is valid
      if [[ " ${genes[*]} " == *" $gene_upper "* ]]; then
        eval "${gene_upper}=true"
      else
        echo "Error: '$gene_name' is not a valid gene name. Valid gene names are: GAG, POL, ENV, VPR, VIF, VPU, NEF." >&2
        exit 1
      fi
  done
fi

# Print chosen genes
listed_genes=()
for gene in "${genes[@]}"; do
  if [ "${!gene}" = true ]; then
    listed_genes+=(" $gene")
  fi
done

if [ "${#listed_genes[@]}" -gt 0 ]; then
  echo -e "${BLUE}Outputting files for chosen options: ${listed_genes[@]} ${END}"
else
  echo "No genes selected for codon correction. Quitting."
  exit 1
fi


### Extract sample and reference sequences for virulign processing
# Extract the first sequence in the file - for MinCov files there will be a second GapsFilled file
Python_Extract_Sample="
from Bio import SeqIO
import sys

input_sequences = SeqIO.parse(sys.argv[1], 'fasta')
sample_seq_file = 'temp_SampleSequence.fasta' 
sample_seq = next(input_sequences)
SeqIO.write([sample_seq], sample_seq_file, 'fasta')
"

"$python2" -c "$Python_Extract_Sample" "$SampleFile" || { echo "CodonCorrection.sh was unable \
to extract the sample sequence. Quitting." >&2; exit 1; }

SampleSequence="temp_SampleSequence.fasta"
SampleSequenceNoQuery="temp_SampleSequence_QueryRemoved.fasta"

# Replace '?' with 'N' within the sequence for BLAST (unable to blast to ? characters)
sed '/^>/!s/?/N/g' "$SampleSequence" > "$SampleSequenceNoQuery"

# BLASTn the target sequence to the possible references
BLASTn_Database="$InitDir/BLASTnDB"
blastn -query "$SampleSequenceNoQuery" -db "$BLASTn_Database"/CC_whole_genome -outfmt \
"10 qseqid qseq evalue sseqid pident qlen qstart qend sstart send frames" -out blastn.out -max_target_seqs 1 \
|| { echo "Failed to execute BLASTn command. Quitting." >&2; exit 1; }
# Change to do top 3-5 then evaluate

blast_output="blastn.out"

# Find the sseqid with the highest pident
# Is this an accurate way to find the closest match? Is evalue better? Also doesn't account for fragment length
RefSequenceName_shell=$(awk -F',' 'NR > 1 {print $5, $0}' "$blast_output" | sort -rn | head -n 1 | awk -F',' '{print $4}')

echo -e "${GREEN}Closest annotatable reference:" "$RefSequenceName_shell${END}"

# Extract the sequence from the reference file
Python_Extract_Reference="
from Bio import SeqIO
import sys

Ref_Name = sys.argv[1]
Ref_File = sys.argv[2]
output_ref = 'temp_Reference_Genome.fasta'

found = False
with open(output_ref, 'w') as file:
  for record in SeqIO.parse(Ref_File, 'fasta'):
    if record.id == Ref_Name:
      SeqIO.write(record, file, 'fasta')
      found = True
      break

if not found:
  print('Reference sequence {} was not found win the provided reference file {}.'.format(Ref_Name, Ref_File))
"

"$python2" -c "$Python_Extract_Reference" "$RefSequenceName_shell" "$ReferenceFile" || { echo "CodonCorrection.sh was unable \
to extract the reference sequence. Quitting." >&2; exit 1; }

# Assign the individual reference file - consists of the full genome sample
ReferenceGenome=temp_Reference_Genome.fasta

###################################################################################

### Align reference to sample
# Creates temp_ files from within the AlignMoreSeqs script that will be overwritten if batch processing 
if [[ "$SingleSequence" == 'true' ]]; then
  # Carry out a mafft pairwise alignment due to no missing coverage in SID_remap_ref files
  # mafft can use the "$mafft" variable if used within shiver
  # options??
  mafft --add "$ReferenceGenome" "$SampleFile" > "$SequenceName_shell"_AlignedSeqs.fasta
else
  # Capable of handling sequence with unknown coverage in parts. Requires two sequences in the starting sample .fasta
  "$python2" "$AlignMoreSeqsTool" "$ReferenceGenome" "$SampleFile" > "$SequenceName_shell"_AlignedSeqs.fasta
fi

ReferenceAlignment="$SequenceName_shell"_AlignedSeqs.fasta

# python also needs splitting into functions
# Add file overwrite checks
# The [...]_Reference_Genes_Extracted.fasta is named after the sample it was generated for, but the
  # individual genes are named after the actual reference. Change this to make it clearer or keep as is?
  # Reason I did it this way was I think it helps to match up the sample genes and the reference used, but
  # once making files for individual genes it seems clearer to separate by actual ID as there's more of them 
  # (also they're intended as temp files)

# Call python script to extract individual genes from the reference and sample
# Set up shiver location properly
CC_Extract_Genes="$HOME/shiver/tools/CC_Extract_Genes.py"

# Check references are divisible by three - i.e. readable as codons

Python_Extract_Genes= # Now removed - clear all references
"$python2" "$CC_Extract_Genes" "$ReferenceAlignment" "$RefSequenceName_shell" "$GeneCoordInfo" "$GAG" "$POL" "$VIF" "$VPR" "$VPU" "$ENV" "$NEF" "$SingleSequence"

###################################################################################

# Determine VIRULIGN options
Nucleotides=false
AminoAcids=false
Mutations=false
MutationTable=false

VirulignOptionsHelp="CodonCorrection.sh can run VIRULIGN with four distinct output formats. To specify \
which, pass up to four of the following options to CodonCorrection.sh enclosed by double quotes and separated by commas.
See the VIRULIGN tutorial for more information on these options.
(1) Nucleotides: An aligned sequence of corrected nucleotides.
(2) AminoAcids: An aligned sequence of corrected amino acids.
(3) Mutations: A list of amino acids changes compared to the reference sequence, including frameshifts.
(4) MutationTable: A CSV file where each mutation present at a specific position is given as a separate column in Boolean representation.
"

# Select VIRULIGN options
VirulignOutput=("Nucleotides" "AminoAcids" "Mutations" "MutationTable")
if [[ "$VirulignOptions" == "all" ]]; then 
  Nucleotides=true
  AminoAcids=true
  Mutations=true
  MutationTable=true
else
  # Select alignment output based on VirulignOptions args
  IFS=',' read -ra listed_options <<< "$VirulignOptions"
  for option in "${listed_options[@]}"; do
    # Could convert to all upper/lowercase here (and change variable names) for conistency
    if [[ " ${VirulignOutput[*]} " == *" $option "* ]]; then # Are * necessary for my data? For this and genes
      eval "${option}=true"
    else
      echo "Error: '$option' is not a valid option. Valid options are Nucleotides, AminoAcids, Mutations, MutationTable. Quitting."
      echo "$VirulignOptionsHelp"
      exit 1
    fi
  done
fi

# Print chosen virulign options
listed_options=()
for option in "${VirulignOutput[@]}"; do
  if [ "${!option}" = true ]; then
    listed_options+=(" $option")
  fi
done
if [ "${#listed_options[@]}" -gt 0 ]; then
  echo -e "${BLUE}Outputting files for chosen options: ${listed_options[@]} ${END}"
else
  echo -e "${RED}Problem with chosen virulign options. Quitting.${END}" >&2
  exit 1
fi

# Override for Nucleotides - in order to correctly analyse length this is always set to true
# TODO should be incorporated more elegantly into the code (currently listed options wont print it, add documentation for being true, 
# adjust options - cannot be removed entirely in the event you only want Nucl as argument cannot be empty)
Nucleotides=true

# Initialise variables
export_alphabet=""
export_kind=""
FileAppend=""

# Define a function to call virulign
function run_virulign {
  export_alphabet_func="$1"
  export_kind_func="$2"
  FileAppend_func="$3"
  # Setting a max value for frameshifts arbitrarily due to frameshift errors occuring in highly variable genes
  # could be adjusted or change based on which gene is being analysed - eg ENV
  MaxFrameshifts=100

  # Find sample and reference files
  for gene in "${genes[@]}"; do
    if [ "${!gene}" == true ]; then
      SampleGene=$(find . -type f -name "temp_${SequenceName_shell}_${gene}_only.fasta")
      ReferenceGene=$(find . -type f -name "temp_${RefSequenceName_shell}_${gene}_Reference_only.fasta")

      if [ -n "$SampleGene" ]; then
        # Run VIRULIGN
        v_output=$( { $VirulignLocation $ReferenceGene $SampleGene --exportKind $export_kind_func --exportAlphabet $export_alphabet_func --exportReferenceSequence yes \
          --exportWithInsertions yes --maxFrameShifts $MaxFrameshifts; } 2>&1 > HIV_${gene}$FileAppend_func )
        echo -e "${GREY}$v_output${END}"

        # Check for errors
        ErrorFile="FailedAlignment.txt"
        touch "$ErrorFile"

        if echo "$v_output" | grep -qi 'error'; then
          error_line=$(echo "$v_output" | grep -i 'error')
          if ! grep -qF "$error_line" "$ErrorFile"; then
            echo "$error_line" >> "$ErrorFile"
          fi
        fi
      else
        echo "Skipping virulign analysis of ${gene} due to missing coverage."
      fi
    fi
  done
}

function call_virulign {
  ### Define virulign command options based on chosen output options
  for option in "${VirulignOutput[@]}"; do
    if [ "${!option}" == "true" ]; then
      # Set variables for each virulign output option
      case "$option" in
        "Nucleotides")
          export_alphabet="Nucleotides"
          export_kind="GlobalAlignment"
          FileAppend="_Nucl_corrected.fasta"
          ;;
        "AminoAcids")
          export_alphabet="AminoAcids"
          export_kind="GlobalAlignment"
          FileAppend="_AminoAcids_corrected.fasta"
          ;;
        "Mutations")
          export_alphabet="Nucleotides"
          export_kind="Mutations"
          FileAppend="_Mutations.csv"
          ;;
        "MutationTable")
          export_alphabet="AminoAcids"
          export_kind="MutationTable"
          FileAppend="_MutationTable.csv"
          ;;
      esac
  
      if [[ -z $export_alphabet || -z $export_kind || -z $FileAppend ]]; then
          echo -e "${RED}Error: One or more options are not set. Quitting.${END}"
          exit 1
      fi

     # Run virulign
     run_virulign "$export_alphabet" "$export_kind" "$FileAppend" || { echo "Problem running virulign. Quitting." >&2; exit 1; }

    fi
  done

  # Check options are set correctly # Possibly no longer useful 
  if [[ -z $export_alphabet || -z $export_kind || -z $FileAppend ]]; then
      echo -e "${RED}Error: One or more options are not set. Quitting.${END}"
      exit 1
  fi
}

# Call virulign with chosen options
call_virulign || { echo "Problem calling virulign. Quitting." >&2; exit 1; } # More descriptive vs run function


# Check for single sequences in output (indicates failed as they should be paired)
if [[ "$Nucleotides" == "true" ]]; then
  for gene in "${genes[@]}"; do
    if [ -f  "HIV_${gene}_Nucl_corrected.fasta" ]; then
      sequence_num=$(grep ">" "HIV_${gene}_Nucl_corrected.fasta" | wc -l)
      if [ "$sequence_num" -lt 2 ]; then
        # Set gene value to false so it's not used for subsequent analysis
        eval "$gene=false"
        echo -e "${YELLOW}Due to failed alignment ${gene} will be omitted from further processing.${END}"
        # If not already in the failed alignments error file, then note the error
        if ! grep -qE "^${SequenceName_shell}_${gene}" "$ErrorFile"; then
          echo "${SequenceName_shell}_${gene} did not produce an alignment due to an unspecified error." >> "$ErrorFile"
        fi
      fi
    fi
  done
fi
# do for other three outputs

# Function to add length to the reference
  # Adds codon-correct number of N's to the reference file (currently just adds to same temp file, could be a separate file)
  # Reasion for the function: For some samples if the reference is shorter than the sample sequence, VIRULIGN will prematurely truncate the sequence
  # to be the same length as the reference. This can be negated by adding indeterminate 'N' bases to the reference. Ideally VIRULIGN would be adjusted 
  # to not require this (I'm unable to find anything within the present options) - but I don't have the C++ code knowledge to correct this.
function add_length {
  LengthChange=$1
  Gene=$2
  ReferenceGene=$(find . -type f -name "temp_${RefSequenceName_shell}_${gene}_Reference_only.fasta")
  NsToAdd=""

  if [ "$LengthChange" -lt 0 ]; then
    echo -e "${YELLOW}The VIRULIGN-corrected $Gene sequence has a change in length of $LengthChange${END}"
    echo "${SequenceName_shell}_$Gene had a change in length of $LengthChange before correction.">> "$LengthFile"
    # Convert to positive number
    LengthChange=$((LengthChange * -1))
    # Determine length of N string to add
    if ((LengthChange % 3 == 0)); then
      N_Number=$LengthChange
    elif ((LengthChange % 3 == 1)); then
      N_Number=$((LengthChange + 2))
    elif ((LengthChange % 3 == 2)); then
      N_Number=$((LengthChange + 1))
    fi
    # Create N string
    for ((i=0; i<N_Number; i++)); do
      NsToAdd+="N"
    done
    echo "$(cat $ReferenceGene)$NsToAdd" > $ReferenceGene

    # Delete old virulign results - could also just be renamed
    rm -f "HIV_${Gene}_Nucl_corrected.fasta"
    rm -f "HIV_${Gene}_AminoAcids_corrected.fasta"
    rm -f "HIV_${Gene}_Mutations.csv"
    rm -f "HIV_${Gene}_MutationTable.csv"
    
    # Set all other gene variables to false in order to re-run virulign with only the modified gene references
    local GAG=false
    local POL=false
    local VIF=false
    local VPR=false
    local VPU=false
    local ENV=false
    local NEF=false

    if [[ "${genes[*]}" == *"$Gene"* ]]; then
      eval "${Gene}=true"
    fi


    # Call virulign again to analyse using new reference lengths 
    call_virulign || { echo "Problem calling virulign after length correction. Quitting." >&2; return 1; }
  fi

}

# Compare sequence length of the extracted sample to final corrected sequence
if [[ "$Nucleotides" == "true" ]]; then
  LengthFile="FailedLengthSeqs.txt"
  for gene in "${genes[@]}"; do
    if [ "${!gene}" = true ]; then
      SampleGene=$(find . -type f -name "temp_${SequenceName_shell}_${gene}_only.fasta")
      # Determine lengths - ignore gaps and N to determine only 'real' changes in lengths, ie false truncation/extension
      if [ -s "HIV_${gene}_Nucl_corrected.fasta" ] && [ -f "$SampleGene" ]; then
        extracted_seq_length=$(awk -v gene="${gene}" -v seq_name="${SequenceName_shell}" '/^>/{if ($0 ~ seq_name "_" gene) found=1; else found=0; next} found { gsub("-", "", $0); len+=length($0) } END{print len}' "$SampleGene") # add error checks to these if failed
        v_output_sample_seq_length=$(awk -v gene="${gene}" -v seq_name="${SequenceName_shell}" '/^>/{if ($0 ~ seq_name "_" gene) found=1; else found=0; next} found { gsub("-", "", $0); gsub("N", "", $0); len+=length($0) } END{print len}' "HIV_${gene}_Nucl_corrected.fasta")
        seq_difference=$((v_output_sample_seq_length - extracted_seq_length))
        if [ "$seq_difference" -lt 0 ]; then
          add_length "$seq_difference" "${gene}" || { echo "Problem adding length to reference files. Quitting." >&2; exit 1; }
        fi
      fi
    fi
  done
fi

# Error message for truncated/extended sequences
if [[ -s $LengthFile ]]; then
  echo -e "${YELLOW}Virulign output files for truncated sample sequences have been deleted. Virulign has been run again with the reference sequence extended for those genes. Ignore the added \
terminal N bases in subsequent analysis - the unmodified reference sequence remains in ${SequenceName_shell}_Reference_GenesExtracted.fasta${END}"
  echo -e "${YELLOW}Check $LengthFile for details of which sequences had their references extended.${END}"
# not sure if this is generated if no length issues; if so, remove it
fi

# improve readibility of printed errors/output

# Check if any alignments failed
if [[ -s "$ErrorFile" ]]; then
  echo -e "${YELLOW}Not all alignments were successful, check $ErrorFile${END}"
else
  rm "$ErrorFile"
fi

# Print frameshift info to file
# TODO if sequence_name already in file return error
# Could be set up to return more meaningful info
if [[ ! -f "Frameshifts.txt"  ]] && [[ "$Mutations" == "true" ]]; then
  echo "Sequence Name, Frameshift quantity, Mutations" > Frameshifts.txt
fi

if [[ -f "Frameshifts.txt" ]] && [[ "$Mutations" == "true" ]]; then
  # Extract frameshift information
  for gene in "${genes[@]}"; do
    if [[ -f "HIV_${gene}_Mutations.csv" ]]; then
      frameshift_info=$(awk -F ',' 'NR==3 {print $1 ", " $4 ", " $7}' "HIV_${gene}_Mutations.csv")
      frameshift_number=$(echo "$frameshift_info" | awk -F', ' '{print $2}')
      # Check if frameshifts are present and that the alignment has not failed
      if [[ "$frameshift_number" != "0" && $(grep -c "${SequenceName_shell}_${gene}" "$ErrorFile") -eq 0 ]]; then
        echo "$frameshift_info" >> Frameshifts.txt
      fi
    fi
  done
elif [[ "$Mutations" == "false" ]]; then
  echo -e "${BLUE}To list frameshifts found in the sample enable the 'Mutations' virulign option.${END}"
fi

# Delete frameshift if no frameshifts are recorded
if [ "$( wc -l <Frameshifts.txt )" -eq 1 ]; then
    rm Frameshifts.txt
    echo "No frameshifts found within successful alignments for $SequenceName_shell"
fi

# Delete temp_ files
# could add a config variable with optional keep temp files, or for some
# rm -f temp_*

# Warn if some genes have missing coverage
if [[ -s "MissingCoverage.txt" ]]; then
  echo -e "${YELLOW}Some genes were not analysed due to missing coverage, check 'MissingCoverage.txt'.${END}"
fi


# After correcting the gene, do we update the full genome consensus?

# Potentially just for testing but check the new sequence BLASTXs
# old code for blastx template below
# blastx -query ~/codon-awareness/Sample/test2_shifted.fasta -db blast_db_pol \
# -outfmt "10 qseqid qseq evalue sseqid pident qlen qstart qend sstart send" \
# -out POL_sequence_blasted.out -max_target_seqs 1

# TODO after: remap to this (i.e. feed into shiver so name files accordingly) verify that 
# the .bam produced is more correct than previously, verify that protein is more correct 
  # could rename the shiver-generated remap file and name my new one the same thing to replace it to preserve shiver code

# clear temp files?? should they not be temp files? some also are overwritten per sample

# Had to set up chmod for first time running. is this a problem for future use in other machines?

# output of virulign: if nucl corrected contains one sequence return as failed

