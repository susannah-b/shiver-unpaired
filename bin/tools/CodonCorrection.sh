#!/usr/bin/env bash

# CodonCorrection.sh is a script designed to call the program 'VIRULIGN' which pairwise aligns viral sequences in a codon-correct manner. Using a set of provided
# references known to be codon-correct and a list of their gene coordinates (see the CodonCorrectionInit.sh), the script will extract the selected genes from the 
# sample and the reference determined to be the closest BLAST match. This is done per individual gene, meaning the reference can be different for each. Each gene
# and the corresponding reference is then fed into VIRULIGN, which will correct the gene sequence by adding 'N' bases to preserve the reading frame in the event of 
# insertions or deletions relative to the reference. The results of this can be viewed in the [Nucl/AminoAcids/Mutations/Mutations]_corrected files depending on 
# the options chosen. If the 'Mutations' option is used, a text file listing frameshifts will be produced for subsequent analysis. 

# CodonCorrectionInit.sh must be run once before running this script

# Explanation of the virulign export options:
# Within this script, there are four different options for the format of the virulign output, which are set in the third argument of the CodonCorrection.sh command.
# Multiple output formats can be generated simultaneously when separated by ',', or use 'all' to generate all four. Each output file will be generated per chosen gene.
  # 1. Nucleotides: will generate a FASTA file of the target nucleotide sequences, with the sample sequence aligned to the reference sequence. This is set to always
  #    be generated as it allows CodonCorrection.sh to check that the sample sequence has not been truncated by VIRULIGN (an issue I observed with some sequences - in 
  #    this event the reference sequence is artificially extended using 'N'.).
  # 2. AminoAcids: will generate a FASTA file of the target amino acid sequences, with the sample sequence aligned to the reference sequence. 
  # 3. Mutations: will output for each gene a list of amino acids changes compared to the reference sequence. This must be enabled in order to generate the file listing
  #    frameshifts.
  # 4. MutationTable: will create a CSV file where each mutation present at a specific position is given as a separate column in Boolean representation. The CSV file is 
  #    annotated according to the numerical position in the protein

# Currently only handles uninterrupted genes GAG POL VIF VPR VPU ENV NEF

set -u
set -o pipefail
set -e

UsageInstructions="Arguments for this script:
(1) A sample file containing the sequence to be corrected as the first sequence;
(2) The filepath for the virulign binaries
(3) The export options for virulign separated by ',' and enclosed by double quotes. One of more of Mutations, MutationTable, Nucleotides, or Amino Acids, \
or use 'all' to generate all four.
(4) The chosen genes to analyse for the sample separated by ',' and enclosed by double quotes. One or more of gag, pol, env, vpr, vif, vpu, or nef, \
or use 'all' to analyse all.
(5) The init directory.
(6) The directory where output files will be created.
(7) OPTIONAL - The number corresponding to how close the BLAST match should be for this run. E.g. the top match is 1, second highest is 2, etc. Default is 1.
"

# Check arguments
ExpectedArgs=6
OptionalArgs=7
if [[ $# -ne "$ExpectedArgs" ]] && [[ "$#" -ne "$OptionalArgs" ]]; then
  echo "$UsageInstructions"
  echo "Incorrect number of arguments provided, expected:" "$ExpectedArgs" "or" "$OptionalArgs" "arguments if selecting which BLAST match number to use." "Quitting." >&2
  exit 1
else
  SampleFile=$1
  VirulignLocation=$2
  VirulignOptions=$3
  GenesToAnalyse=$4
  InitDir=$5
  OutputDir=$6
  if [[ $# -eq 6 ]]; then
    BLASTMatch=1
  elif [[ $# -eq 7 ]]; then
    BLASTMatch=$7
  fi
fi

# Find shiver files
ToolsDir="$(cd "$(dirname "${BASH_SOURCE[0]}" )" && pwd)"
shiver="$(dirname "$ToolsDir")" # shiver file location
source "$shiver"/'bin/config.sh'
PythonFuncs="$shiver"/'bin/tools/CC_Python_funcs.py'
ReferenceFile="$InitDir/CC_References.fasta"
GeneCoordInfo="$InitDir/CC_Coords.fasta"
HXB2File="$shiver"/'data/external/B.FR.83.HXB2_LAI_IIIB_BRU.K03455.fasta'
AlignMoreSeqsTool="$shiver"/'bin/tools/AlignMoreSeqsToPairWithMissingCoverage.py'

# Colours # Can be removed as they don't work in the cluster. For now I've set a variable in the config to disable them
if [[ "$EnableColours" == "true" ]]; then
  RED="\033[91m"
  BLUE="\033[94m"
  YELLOW="\033[93m"
  GREEN="\033[92m"
  GREY="\033[37m"
  END="\033[0m"
else
  RED=""
  BLUE=""
  YELLOW=""
  GREEN=""
  GREY=""
  END=""
fi

# Check init script has been run
if [[ ! -s "$ReferenceFile" ]] || [[ ! -s "$GeneCoordInfo" ]]; then
  echo "Reference file and/or the coordinates file are not found in the init directory. Please run CodonCorrectionInit.sh before running this script. Quitting." >&2
  exit 1
fi

# Check if file exists
if ! [[ -f "$SampleFile" ]]; then
  echo "The provided sample file '$SampleFile' does not exist. Check the file path is correct."
  exit 1
fi

# Check sample is a .fasta file
if [[ "$SampleFile" != *.fasta ]]; then 
  echo "Sample file $SampleFile is not a fasta file. Quitting." >&2
  exit 1
fi

# Count sequence number in sample 
Sample_SeqNumber=$(grep '^>' "$SampleFile" | wc -l | awk '{$1=$1};1')
if [[ "$Sample_SeqNumber" == 0 ]]; then
  echo "Reference file $ReferenceFile contains no sequences. Quitting." >&2
  exit 1
fi

# Check HXB2 file exists within shiver
if ! [[ -f "$HXB2File" ]]; then
  echo "Unable to find the HXB2 file within shiver: $HXB2File. Quitting." >&2
  exit 1
fi

# Determine whether to use mafft pairwise or AlignMoreSeqs
if [[ "$Sample_SeqNumber" == 1 ]]; then
  SingleSequence=true # Will be funneled into mafft pairwise alignment
elif [[ "$Sample_SeqNumber" == 2 ]]; then
  SingleSequence=false # In this case CC.sh expects remap consensus files (=2 sequences) and will use the AlignMoreSeqs tool which can handle missing coverage.
else 
  echo "Unexpected number of sequences in the sample file. Expected 1 or 2. Quitting." >&2
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
{ echo "Unable to change directory to $OutputDir. Quitting." >&2 ; exit 1 ; }

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
    # Could convert to all upper/lowercase here (and change variable names) for consistency
    if [[ " ${VirulignOutput[*]} " == *" $option "* ]]; then
      eval "${option}=true"
    else
      echo "Error: '$option' is not a valid option. Valid options are Nucleotides, AminoAcids, Mutations, MutationTable. Quitting."
      echo "$VirulignOptionsHelp"
      exit 1
    fi
  done
fi

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

# Extracting the sample ID from the sample fasta file
SequenceName_shell=$(head -n 1 "$SampleFile" | sed 's/>//')

###################################################################################

# Print chosen genes
listed_genes=()
for gene in "${genes[@]}"; do
  if [[ "${!gene}" = true ]]; then
    listed_genes+=(" $gene")
  fi
done

if [[ "${#listed_genes[@]}" -gt 0 ]]; then
  echo -e "${BLUE}Outputting files for chosen options: ${listed_genes[@]} ${END}"
else
  echo "No genes selected for codon correction. Quitting."
  exit 1
fi

# Function to align the sample to the chosen reference
  # Both mafft options are set to quiet to avoid obfuscating error messages, but this can be switched off by removing the '---quiet' options below
function align_to_sample {
  AlignmentRef="$1"
  AlignmentAppend="$2"

  if [[ "$AlignmentRef" == "$HXB2File" ]]; then
    OutputAlignment=temp_"$SequenceName_shell$AlignmentAppend"
  else
    OutputAlignment="$SequenceName_shell$AlignmentAppend"
  fi

  if [[ "$SingleSequence" == 'true' ]]; then
    # Carry out a mafft pairwise alignment - this assumes that single sequences do not contains the '?' characters found within shiver output files 
    "$mafft" --quiet --add "$AlignmentRef" "$SampleFile" > "$OutputAlignment"
  else
    # Capable of handling sequence with unknown coverage in parts. Requires two sequences in the starting sample .fasta
    "$python" "$AlignMoreSeqsTool" --x-mafft "$mafft --quiet" "$AlignmentRef" "$SampleFile" > "$OutputAlignment"
  fi

}

# Align sample to HXB2 in preparation for gene extraction
align_to_sample "$HXB2File" "_HXB2Aligned.fasta" || { echo "Problem aligning sample with reference HXB2. Quitting." >&2; exit 1; }
HXB2_Alignment=temp_"$SequenceName_shell$AlignmentAppend"

# Extract approximate sample genes using HXB2 as a reference - only used to BLAST to reference genes before more precise sample extraction
"$python" "$PythonFuncs" ExtractWithHXB2 --AlignmentFile "$HXB2_Alignment" --SingleSeq "$SingleSequence" || { echo "CodonCorrection.sh was unable to extract the sample gene sequences. Quitting." >&2; exit 1; }

# Check if genes have entered MissingCoverage.txt, and if so omit from further processing
if [[ -s "MissingCoverage.txt" ]]; then
  for gene in "${genes[@]}"; do
    if [[ "${!gene}" = true ]]; then
      if grep -q "$SequenceName_shell"_"${gene}" "MissingCoverage.txt"; then
          eval "$gene=false"
      fi
    fi
  done
fi

# BLASTn the sample sequences to the corresponding gene database
function blastn_to_genes() {
  MaxTargets="$BLASTMatch" # Ensures that enough targets are listed - can be set to a static value instead of an arg, or a minimum of N.
  for gene in "${genes[@]}"; do
    if [[ "${!gene}" = true ]]; then
      BLASTn_Database="$InitDir/BLASTnDB_${gene}"
      BLAST_Gene=/"$OutputDir"/temp_preBLAST_Sample_Gene_"${gene}".fasta
      BLAST_Output=blastn_"${gene}".out
      # Run BLASTn
      "$BlastNcommand" -query "$BLAST_Gene" -db "$BLASTn_Database"/"${gene}" -outfmt \
      "10 qseqid qseq evalue sseqid pident qlen qstart qend sstart send bitscore" -out "$BLAST_Output" -max_target_seqs "$MaxTargets" \
      || { echo "Failed to execute BLASTn command. Quitting." >&2; return 1; }
      
      # Assign the nth closest reference name per gene, where n is the specified nth closest match (1 is 1st, 2 is 2nd, etc.)
      ExtractMatch="$(awk -F',' '{print $11, $0}' "$BLAST_Output" | sort -rn | awk -F',' -v rank="$BLASTMatch" 'NR==rank {print $4}' | sed 's/^[[:space:]]*//')"
      ExtractSeqName=${ExtractMatch%????}

      # Assign a reference variable for each gene
      GeneRef="ref_$gene"
      eval "$GeneRef=$ExtractSeqName"
      # List reference to shell
      echo -e "${GREEN}Reference used for ${gene}: ${!GeneRef}${END}"
    fi
  done
  return 0
}

# Call BLASTn function
blastn_to_genes || { echo "Problem BLASTing sample genes to the corresponding gene database. Quitting." >&2; exit 1; }

# Extract the reference sequence from the full references list for each gene
for gene in "${genes[@]}"; do
  if [[ "${!gene}" = true ]]; then
    GeneRef="ref_$gene"
    "$python" "$PythonFuncs" ExtractRefFromFasta --ReferenceName "${!GeneRef}" --GenomeFile "$ReferenceFile" --Gene "$gene" || { echo "CodonCorrection.sh was unable to extract the reference sequence from the list of references. Quitting." >&2; exit 1; }
  fi
done

# Assign the reference genome as a variable for each gene
for gene in "${genes[@]}"; do
  if [[ "${!gene}" = true ]]; then
    ReferenceGenome="ref_genome_$gene"
    eval "$ReferenceGenome=temp_Reference_Genome_$gene.fasta"
  fi
done

# Align the sample genome to the reference genome for each gene and extract the sample and reference genes
  # Sample genes have been extracted with HXB2 already (although perhaps less accurately as it just uses HXB2 to align and extract) 
  # and the reference genes were extracted as part of the init. Because this script is already written I'm keeping it for the updated 
  # gene-by-gene processing, but it could be modified/removed due to those preexisting files.
for gene in "${genes[@]}"; do
  if [[ "${!gene}" = true ]]; then
    ReferenceGenome="ref_genome_$gene"
    GeneRef="ref_$gene"
    align_to_sample "${!ReferenceGenome}" "_AlignedSeqs_$gene.fasta" || { echo "Problem aligning sample with reference genome. Quitting." >&2; exit 1; }
    # Assign the reference alignment file
    ReferenceAlignment="$SequenceName_shell"_AlignedSeqs_$gene.fasta # This variable will only work within the loop

   # Run python script to extract the final sample and reference genes
     # Names of output files could be confusing? Rename if needed
   "$python" "$PythonFuncs" ExtractRefandSample --AlignmentToRef "$ReferenceAlignment" --SingleSeq "$SingleSequence" --ReferenceName "${!GeneRef}" --GeneCoordInfo "$GeneCoordInfo" --Gene "${gene}" || { echo "CodonCorrection.sh was unable to extract the genes for reference and sample. Quitting." >&2; exit 1; }
  fi
done

###################################################################################

# Print chosen virulign options
listed_options=()
for option in "${VirulignOutput[@]}"; do
  if [[ "${!option}" = true ]]; then
    listed_options+=(" $option")
  fi
done
if [[ "${#listed_options[@]}" -gt 0 ]]; then
  echo -e "${BLUE}Outputting files for chosen options: ${listed_options[@]} ${END}"
else
  echo -e "${RED}Problem with chosen virulign options. Quitting.${END}" >&2
  exit 1
fi

# Override for Nucleotides - in order to correctly analyse length this is always set to true
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
  # Setting a max value for frameshifts arbitrarily due to frameshift errors occuring in highly variable genes - although in the gene I tested the same "Frameshift Error"
  # occured within VIRULIGN regardless of this value.
  MaxFrameshifts=50
  GapExtensionPenalty=3.3 # default virulign is 3.3
  GapOpenPenalty=10.0 # default virulign is 10.0

  # Find sample and reference files
  for gene in "${genes[@]}"; do
    if [[ "${!gene}" == true ]]; then
      GeneRef="ref_$gene"
      SampleGene=$(find . -type f -name "temp_${SequenceName_shell}_${gene}_only.fasta")
      ReferenceGene=$(find . -type f -name "temp_${!GeneRef}_${gene}_Reference_only.fasta")

      if [[ -n "$SampleGene" ]]; then
        # Run VIRULIGN
        v_output=$( { $VirulignLocation $ReferenceGene $SampleGene --exportKind $export_kind_func --exportAlphabet $export_alphabet_func --exportReferenceSequence yes \
          --exportWithInsertions yes --gapExtensionPenalty $GapExtensionPenalty --gapOpenPenalty $GapOpenPenalty --maxFrameShifts $MaxFrameshifts; } 2>&1 > HIV_${gene}$FileAppend_func )
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

# Function to call the run_virulign function with the chosen VIRULIGN options
function call_virulign {
  # Define virulign command options based on chosen output options
  for option in "${VirulignOutput[@]}"; do
    if [[ "${!option}" == "true" ]]; then
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
}

# Call virulign with chosen options
call_virulign || { echo "Problem calling virulign with the selected options. Quitting." >&2; exit 1; }


# Check for single sequences in output (indicates failed as they should be paired)
for gene in "${genes[@]}"; do
  # Check Nucleotides aligned successfully
  if [[ "$Nucleotides" == "true" ]]; then
    if [[ -f  "HIV_${gene}_Nucl_corrected.fasta" ]]; then
      sequence_num_nucl=$(grep ">" "HIV_${gene}_Nucl_corrected.fasta" | wc -l)
      if [[ "$sequence_num_nucl" -lt 2 ]]; then
        # Set gene value to false so it's not used for subsequent analysis
        eval "$gene=false"
        echo -e "${YELLOW}Due to failed nucleotide alignment ${gene} will be omitted from further processing.${END}"
        # If not already in the failed alignments error file, then note the error
        if ! grep -qE "^${SequenceName_shell}_${gene}" "$ErrorFile"; then
          echo "${SequenceName_shell}_${gene} failed to produce an alignment." >> "$ErrorFile"
        fi
        continue
      fi
    fi
  fi
  # Check amino acids aligned successfully
  if [[ "$AminoAcids" == "true" ]]; then
    if [[ -f  "HIV_${gene}_AminoAcids_corrected.fasta" ]]; then
      sequence_num_AA=$(grep ">" "HIV_${gene}_AminoAcids_corrected.fasta" | wc -l)
      if [[ "$sequence_num_AA" -lt 2 ]]; then
        # Set gene value to false so it's not used for subsequent analysis
        eval "$gene=false"
        echo -e "${YELLOW}Due to failed amino acid alignment ${gene} will be omitted from further processing.${END}"
        # If not already in the failed alignments error file, then note the error
        if ! grep -qE "^${SequenceName_shell}_${gene}" "$ErrorFile"; then
          echo "${SequenceName_shell}_${gene} failed to produce an alignment." >> "$ErrorFile"
        fi
        continue
      fi
    fi
  fi
  # Check Mutations file analysed both sequences 
  if [[ "$Mutations" == "true" ]]; then
    # Extract frameshift information
    if [[ -f "HIV_${gene}_Mutations.csv" ]]; then
      MutationCheck=$(awk -F ',' 'NR==3 {print $2}' "HIV_${gene}_Mutations.csv")
      # Check if the alignment failed
      if [[ "$MutationCheck" == "Failure" ]]; then
        # Set gene value to false so it's not used for subsequent analysis
        eval "$gene=false"
        echo -e "${YELLOW}Due to failed nucleotide alignment ${gene} will be omitted from further processing.${END}"
        # If not already in the failed alignments error file, then note the error
        if ! grep -qE "^${SequenceName_shell}_${gene}" "$ErrorFile"; then
          echo "${SequenceName_shell}_${gene} failed to produce an alignment." >> "$ErrorFile"
        fi
        continue
      fi
    fi
  fi
  # Check Mutation table contains both sequences
  if [[ "$MutationTable" == "true" ]]; then
    # Count lines
    if [[ -f "HIV_${gene}_MutationTable.csv" ]]; then
      mt_line_count=$(cat "HIV_${gene}_MutationTable.csv" | wc -l)
      if [[ "$mt_line_count" -eq 2 ]]; then
        # Set gene value to false so it's not used for subsequent analysis
        eval "$gene=false"
        echo -e "${YELLOW}Due to failed amino acid alignment ${gene} will be omitted from further processing.${END}"
        # If not already in the failed alignments error file, then note the error
        if ! grep -qE "^${SequenceName_shell}_${gene}" "$ErrorFile"; then
          echo "${SequenceName_shell}_${gene} failed to produce an alignment." >> "$ErrorFile"
        fi
        continue  
      fi
    fi
  fi
done

# Function to add length to the reference
  # Adds codon-correct number of N's to the temp reference file. For some samples if the reference is shorter than the sample sequence, VIRULIGN will prematurely 
  # truncate the sequence to be the same length as the reference. This can be negated by adding generic 'N' bases to the reference. Ideally VIRULIGN would 
  # be adjusted to not require this (I'm unable to find anything within the present options) - but I don't have the C++ knowledge to correct this.
function add_length {
  LengthChange=$1
  Gene=$2
  ReferenceGene=$(find . -type f -name "temp_${!GeneRef}_${gene}_Reference_only.fasta")
  NsToAdd=""

  # Extract the sample sequence from the first round of virulign correction
  "$python" "$PythonFuncs" ExtractSequence --OutputFile "temp_${SequenceName_shell}_${Gene}_Corrected.txt" --InputFile "HIV_${Gene}_Nucl_corrected.fasta" --SequenceNumber "2" || { echo "Problem extracting the corrected gene sequence. Quitting." >&2; return 1; }
  # Determine if any bases are missing from the start of the virulign output by aligning the extracted and corrected gene
  SampleAlignment="temp_${SequenceName_shell}_${Gene}_SampleAlignment.fasta"
  "$mafft" --quiet --add "temp_${SequenceName_shell}_${Gene}_only.fasta" "temp_${SequenceName_shell}_${Gene}_Corrected.txt" > "$SampleAlignment"

  # Count the gaps at the start of the corrected sequence (i.e. the shift in gene positions of the sample sequence)
  PositionShift=0
  PositionShift=$(awk 'NR==2 {if(match($0, /^[-]+/)) print RLENGTH; else print 0}' "$SampleAlignment")
  if [[ "$LengthChange" -lt 0 ]]; then
    echo -e "${YELLOW}The VIRULIGN-corrected $Gene sequence has a change in length of $LengthChange${END}"
    echo "${SequenceName_shell}_$Gene had a change in length of $LengthChange before correction.">> "$LengthFile"
    # Convert to positive number (length change is always negative for truncation issues)
    LengthChange=$((LengthChange * -1))
    # Determine change in length at the end of the sequence
    EndLength=$((LengthChange - PositionShift))
    if [[ "$EndLength" != 0 ]]; then
      # Determine length of N string to add
      if ((EndLength % 3 == 0)); then
        N_Number=$EndLength
      elif ((EndLength % 3 == 1)); then
        N_Number=$((EndLength + 2))
      elif ((EndLength % 3 == 2)); then
        N_Number=$((EndLength + 1))
      fi
      # Create N string
      for ((i=0; i<N_Number; i++)); do
        NsToAdd+="N"
      done
      echo "$(cat $ReferenceGene)$NsToAdd" > "$ReferenceGene"
    fi

    # Add length to the start of the sequence
    if [[ "$PositionShift" != 0 ]]; then
      NsToAdd=""
      for ((i=0; i<PositionShift; i++)); do
        NsToAdd+="N"
      done
      awk 'NR==2 {print "'"$NsToAdd"'" $0; next} 1' "$ReferenceGene" > temp_file && mv temp_file "$ReferenceGene"
    fi

    # OPTION 1: Delete old results
    # Delete old virulign results - For testing purposes files are renamed as below, can revert to deletion
    # rm -f ./HIV_${Gene}_Nucl_corrected.fasta
    # rm -f ./HIV_${Gene}_AminoAcids_corrected.fasta
    # rm -f ./HIV_${Gene}_Mutations.csv
    # rm -f ./HIV_${Gene}_MutationTable.csv

    # OPTION 2: Rename virulign results
    mv -f ./HIV_${Gene}_Nucl_corrected.fasta ./HIV_${Gene}_Nucl_corrected_OLD.fasta
    mv -f ./HIV_${Gene}_AminoAcids_corrected.fasta ./HIV_${Gene}_AminoAcids_corrected_OLD.fasta
    mv -f ./HIV_${Gene}_Mutations.csv ./HIV_${Gene}_Mutations_OLD.csv
    mv -f ./HIV_${Gene}_MutationTable.csv ./HIV_${Gene}_MutationTable_OLD.csv
    
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
  LengthFile="CorrectedLengthSeqs.txt"
  for gene in "${genes[@]}"; do
    GeneRef="ref_$gene"
    if [[ "${!gene}" = true ]]; then
      SampleGene=$(find . -type f -name "temp_${SequenceName_shell}_${gene}_only.fasta")
      # Determine lengths - ignore gaps and added Ns to determine only 'real' changes in lengths, ie false truncation/extension
      if [[ -s "HIV_${gene}_Nucl_corrected.fasta" ]] && [[ -f "$SampleGene" ]]; then
        extracted_seq_length=$(awk -v gene="${gene}" -v seq_name="${SequenceName_shell}" '/^>/{if ($0 ~ seq_name "_" gene) found=1; else found=0; \
        next} found { gsub("-", "", $0); len+=length($0) } END{print len}' "$SampleGene") # add error checks to these if failed
        # Find the number of pre-existing N's in the sample gene
        ExistingNs=$(tail -n +2 "$SampleGene" | tr -cd 'Nn' | wc -c | awk '{print $1+0}')
        v_output_sample_seq_length=$(awk -v gene="${gene}" -v seq_name="${SequenceName_shell}" '/^>/{if ($0 ~ seq_name "_" gene) found=1; else found=0; \
        next} found { gsub("-", "", $0); gsub("N", "", $0); len+=length($0) } END{print len}' "HIV_${gene}_Nucl_corrected.fasta")
        v_output_sample_seq_length=$((v_output_sample_seq_length + ExistingNs))
        # Calculate difference in sequence length and add length if needed
        seq_difference=$((v_output_sample_seq_length - extracted_seq_length))
        if [[ "$seq_difference" -lt 0 ]]; then
          add_length "$seq_difference" "${gene}" || { echo "Problem adding length to reference files. Quitting." >&2; exit 1; }
        fi
      fi
    fi
  done
fi

# Error message for truncated/extended sequences
if [[ -s $LengthFile ]]; then
  echo -e "${YELLOW}Virulign output files for truncated sample sequences have been deleted. Virulign has been run again with the reference sequence extended for those genes. Ignore the added \
terminal N bases in subsequent analysis - the unmodified reference sequence remains in ${SequenceName_shell}_Reference_GenesExtracted.fasta.
As a result of this some positional/mutation information for the reference sequence may be incorrect due to artifical bases being added - eg X mutations in the reference.${END}"
  echo -e "${YELLOW}Check $LengthFile for details of which sequences had their references extended.${END}"
fi

# Check if any alignments failed
if [[ -s "$ErrorFile" ]]; then
  echo -e "${YELLOW}Not all alignments were successful, check $ErrorFile${END}"
else
  rm "$ErrorFile"
fi

# Print frameshift info to file
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
      if [[ "$frameshift_number" != "0" && (! -f "$ErrorFile" || $(grep -c "${SequenceName_shell}_${gene}" "$ErrorFile") -eq 0) ]]; then
        echo "$frameshift_info" >> Frameshifts.txt
      fi
    fi
  done
elif [[ "$Mutations" == "false" ]]; then
  echo -e "${BLUE}To list frameshifts found in the sample enable the 'Mutations' virulign option.${END}"
fi

# Delete frameshift if no frameshifts are recorded
if [[ -f "Frameshifts.txt"  ]] && [[ "$Mutations" == "true" ]]; then
  if [[ "$( wc -l <Frameshifts.txt )" -eq 1 ]]; then
      rm Frameshifts.txt
      echo "No frameshifts found within successful alignments for $SequenceName_shell"
  fi
fi

# Identify indels
if [[ -f "Frameshifts.txt" ]]; then
  OutputCSV="Frameshifts.csv"
  for gene in "${genes[@]}"; do
  GeneRef="ref_$gene"
  SampleGene=$(find . -type f -name "temp_${SequenceName_shell}_${gene}_only.fasta")
    if [[ "${!gene}" = true ]]; then
      "$python" "$PythonFuncs" CategoriseIndels --InputFile "HIV_${gene}_Nucl_corrected.fasta" --OutputCSV "$OutputCSV" --Gene "${gene}" \
      --SequenceName "${SequenceName_shell}" --ReferenceName "${!GeneRef}" --GeneCoordInfo "$GeneCoordInfo" --GenomeFile "$ReferenceFile" \
      --AlignmentFile "$HXB2_Alignment" --ExtractedGeneFile "$SampleGene" || { echo "Problem finding indel positions. Quitting." >&2; exit 1; }
    fi
  done
fi

# Warn if some genes have missing coverage
if [[ -s "MissingCoverage.txt" ]]; then
  echo -e "${YELLOW}Some genes were not analysed due to missing coverage, check 'MissingCoverage.txt'.${END}"
fi

# Delete temporary files
# Check that the directory is correct before deleting
if [[ "$DeleteTemp" == "true" ]]; then
  CurrentDir=$(pwd)
  if [[ "$CurrentDir" = "$OutputDir" ]]; then
    rm ./temp_*
    echo "Deleting temporary files. Set DeleteTemp to false in config if choosing to keep these files."
  else
    echo "The current directory is not the expected Output Directory '$OutputDir'. Aborting deletion."
  fi
else
  echo "Temporary files have not been deleted. To do so, set DeleteTemp to true in the config."
fi
