#!/usr/bin/env bash

# CodonCorrectionInit.sh must be run once before using this script

set -u
set -o pipefail
set -e

UsageInstructions="Arguments for this script:
(1) A sample file containing the sequence to be corrected as the first sequence;
(2) A reference file containing all codon-correct/annotatable references;
(3) A file of gene coordinates - NB this requires specific formatting/order (see note in code)
(4) The filepath for the virulign binaries
(5) The export options for virulign separated by ',' and enclosed by double quotes. One of more of Mutations, MutationTable, \
Nucleotides, or Amino Acids, or use 'all' to generate all four.
(6) The chosen genes to analyse for the sample separated by ',' and enclosed by double quotes. \
One or more of gag, pol, env, vpr, vif, vpu, or nef, or use 'all' to analyse all.
"

# Write help text to explain the different export options

### Check arguments
ExpectedArgs=6
if [ $# -ne "$ExpectedArgs" ]; then
  echo "$UsageInstructions"
  echo "Incorrect number of arguments provided, expected:" "$ExpectedArgs". "Quitting." >&2
  exit 1
else
  SampleFile=$1
  ReferenceFile=$2
  GeneCoordInfo=$3
  VirulignLocation=$4
  VirulignOptions=$5 # VirulignOptions and GenesToAnalyse aren't actually used in the rest of the code so far, only $5 and $6
  GenesToAnalyse=$6
fi

# Check init has been done
BLASTnDB_Path=~/shiver/CodonCorrectionReferences/BLASTnDB

if [[ ! -d "$BLASTnDB_Path" ]]; then
  echo "No BLASTn database found. Please run CodonCorrectionInit.sh before running this script. Quitting." >&2
  exit 1
fi

# Check each arg is valid and as expected
  # sample file should contain at least 1 sequence and be a .fasta
  # reference file should be a .fasta
  # check gene coordinates is a .txt with a header of data in the right order
  # test virulign

# Check VIRULIGN exists
if [[ ! -f "$VirulignLocation" ]]; then
  echo "The specified VIRULIGN location '$VirulignLocation' does not exist. Check file path is correct. Quitting." >&2
  exit 1
fi

# IMPORTANT: Gene Coords must be set up in the correct format. My scripts use gene fragments to determine gene coordinates
# from complete genomes and make a specific .txt, however I haven't integrated this into the CodonCorrection script. 
# If just providing a set of references with this tool, then perhaps parts could be simplified to reflect this

### Other variables
# I think shiver assumes it's in ~/shiver so this would also work..? 
AlignMoreSeqsTool=~/shiver/tools/AlignMoreSeqsToPairWithMissingCoverage.py
# Extracting the sample ID from the fasta
SequenceName_shell=$(head -n 1 "$SampleFile" | sed 's/>//')

# Temp for testing
cd ~/shiver
# Create directory to store results
if [ ! -e "CodonCorrection" ]; then
  mkdir CodonCorrection
else
  echo "Using CodonCorrection folder - testing only" # temporary
  # should be proper error message/handling to not overwrite
  # addend SID? or other method to handle identical files for batch processing
  # echo "CodonCorrection file already exists, empty and delete to avoid overwriting."
  # exit 1
  # Could use separate user-defined identifier in command eg shiver
fi

# Not necessary for function but keeps results contained - can be changed
cd CodonCorrection

# Check config and other shiver scripts for variables to work from if needed

###################################################################################

# Check if file exists
if ! [ -f "$SampleFile" ]; then
  echo "The provided sample file does not exist. Check the file path is correct."
  exit 1
# else check .fasta
# also check other files
# not necessarily for this one but check some have right amount of sequences etc
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
if [[ "$6" == "all" ]]; then
  GAG=true
  POL=true
  VIF=true
  VPR=true
  VPU=true
  ENV=true
  NEF=true
else
  # Switch on the corresponding gene variable for each listed gene
  IFS=',' read -ra listed_genes <<< "$6"
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
  echo "Outputting files for chosen options: ${listed_genes[@]}"
else
  echo "No genes selected for codon correction. Quitting."
  exit 1
fi

# Temporary set variable for remap testing - incorporate as argument/auto-detect later
# variable name tbd
# set to false because mafft pairwise align needs to be set up, AlignTool requires 2 seqs
Remap=false


### Extract sample and reference sequences for virulign processing
# Extract the first sequence in the file - for MinCov files there will be a second GapsFilled file (ignore this)
Python_Extract_Sample="
from Bio import SeqIO
import sys

input_sequences = SeqIO.parse(sys.argv[1], 'fasta')
sample_seq_file = 'temp_SampleSequence.fasta' 
sample_seq = next(input_sequences)
SeqIO.write([sample_seq], sample_seq_file, 'fasta')
"

python2 -c "$Python_Extract_Sample" "$SampleFile"

SampleSequence="temp_SampleSequence.fasta"

# BLASTn the target sequence to the possible references
blastn -query "$SampleSequence" -db ~/shiver/CodonCorrectionReferences/BLASTnDB/CC_whole_genome -outfmt \
"10 qseqid qseq evalue sseqid pident qlen qstart qend sstart send frames" -out blastn.out -max_target_seqs 1 
# Change to do top 3-5 then evaluate

blast_output="blastn.out"

# Find the sseqid with the highest pident
# Is this an accurate way to find the closest match? Is evalue better? Also doesn't account for fragment length
RefSequenceName_shell=$(awk -F',' 'NR > 1 {print $5, $0}' "$blast_output" | sort -rn | head -n 1 | awk -F',' '{print $4}')

echo "Closest annotatable reference:" "$RefSequenceName_shell"

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

python2 -c "$Python_Extract_Reference" "$RefSequenceName_shell" "$ReferenceFile"

# Assign the individual reference file - consists of the full genome sample
ReferenceGenome=temp_Reference_Genome.fasta

###################################################################################

### Align reference to sample
# Creates temp_ files which within the AlignMoreSeqs script that will be overwritten if batch processing 
# Can just be a simple pairwise mafft align for remap files (use Remap variable to determine)
python2 "$AlignMoreSeqsTool" "$ReferenceGenome" "$SampleFile" > "$SequenceName_shell"_AlignedSeqs.fasta

ReferenceAlignment="$SequenceName_shell"_AlignedSeqs.fasta

# Definitely better ways to use python within bash but this is easier for testing currently
  # Will create a separate .py later in development
# Also needs splitting into functions
# Add file overwrite checks
# The [...]_Reference_Genes_Extracted.fasta is named after the sample it was generated for, but the
  # individual genes are named after the actual reference. Change this to make it clearer or keep as is?
  # Reason I did it this way was I think it helps to match up the sample genes and the reference used, but
  # once making files for individual genes it seems clearer to separate by actual ID as there's more of them 
  # (also they're intended as temp files)

Python_Extract_Genes="""
import sys
from Bio import AlignIO
from Bio import SeqIO

# Determine position of reference sequence in alignment file
# add check for RefSequenceNumber is not 0 later
Remap_check = sys.argv[11]
RefSequenceNumber = 0
if Remap_check == 'true':
  RefSequenceNumber = 1
else:
  RefSequenceNumber = 2

# Assign variables
ReferenceAlignment = AlignIO.read(sys.argv[1], 'fasta')
SequenceName = ReferenceAlignment[0].id
RefSequenceLength = len(ReferenceAlignment[RefSequenceNumber])
RefSequenceName = sys.argv[2]
GeneCoordInfo = sys.argv[3]
GAG_check = sys.argv[4]
POL_check = sys.argv[5]
VIF_check = sys.argv[6]
VPR_check = sys.argv[7]
VPU_check = sys.argv[8]
ENV_check = sys.argv[9]
NEF_check = sys.argv[10]
output_file_name = '{}_GenesExtracted.fasta'.format(SequenceName)
ref_output_file_name = '{}_Reference_GenesExtracted.fasta'.format(SequenceName)



# The following python code contains lots of gap variables etc, most of this actually isn't
# needed so I removed the parts that adjust gene length based on this, but the rest of the
# code still remains in the script. Can be later removed and some parts (eg enumerate) consensus_num_gaps
# can be simplified (similar to the reference gene extraction that uses sequence length)

# Extract gene coordinates from the reference file
gene_loci = {}
with open(GeneCoordInfo, 'r') as file:
  for line in file:
    gene_info = line.strip().split(',')
    if gene_info[0] == RefSequenceName:
      gene_loci = {
        # Gene: (gene_start, gene_end)
        'GAG': (int(gene_info[1]), int(gene_info[2])),
        'POL': (int(gene_info[3]), int(gene_info[4])),
        'VIF': (int(gene_info[5]), int(gene_info[6])),
        'VPR': (int(gene_info[7]), int(gene_info[8])),
        'VPU': (int(gene_info[9]), int(gene_info[10])),
        'ENV': (int(gene_info[11]), int(gene_info[12])),
        'NEF': (int(gene_info[13]), int(gene_info[14]))
      }
      break

# Gene Extraction for reference sequence
if gene_loci:
  with open(ref_output_file_name, 'w') as output_file:
    for gene, (gene_start, gene_end) in gene_loci.items():
      check_gene = gene + '_check'
      if locals().get(check_gene) == 'true':
        if gene_start is not None and gene_end is not None:
          ref_pos = 0
          ref_start = 0
          ref_end = 0
          for i in xrange(RefSequenceLength):
            if ReferenceAlignment[RefSequenceNumber][i] != '-':
              ref_pos += 1
            if ref_pos == gene_start:
              ref_start = i
            if ref_pos == gene_end:
              ref_end = i
              break
          if ref_start is not None and ref_end is not None:
            reference_sequence = ReferenceAlignment[RefSequenceNumber].seq[ref_start:ref_end + 1].ungap('-')
            output_file.write('>' + ReferenceAlignment[RefSequenceNumber].id + '_' + gene + '\n' + str(reference_sequence) + '\n')

          # testing
          print('ref_start:', ref_start)
          print('ref_end:', ref_end)
          print('gene_start:', gene_start)
          print('gene_end:', gene_end)
          print('ref_pos:', ref_pos)


# Gene extraction for sample sequence
# Select which genes to extract
if gene_loci:
  with open(output_file_name, 'w') as output_file:
    for gene, (gene_start, gene_end) in gene_loci.items():
      check_gene = gene + '_check'
      if locals().get(check_gene) == 'true':
        # Iterate over the sequence to determine gene coordinates of the sample sequence
        consensus_start = None
        consensus_end = None
        sequence_pos = 0
        indel_difference = 0
        for i in xrange(RefSequenceLength):
          if ReferenceAlignment[RefSequenceNumber][i] != '-':
            sequence_pos += 1
          if sequence_pos == gene_start:
            consensus_start = i
          if sequence_pos == gene_end: 
            consensus_end = i
            break
        indel_difference =  (consensus_end + 1 - consensus_start) - (gene_end + 1 - gene_start)
        if indel_difference != 0:
          print ('{} Indel(s) detected in {} for {}.'.format(indel_difference, SequenceName, gene))
        
        #testing
        print('consensus_start:', consensus_start)
        print('consensus_end:', consensus_end)
        print('sequence_pos:', sequence_pos)
        print('indel_difference:', indel_difference)

        # Write the gene sequence to file unless it contains '?'
        if consensus_start is not None and consensus_end is not None:
          gene_sequence = ReferenceAlignment[0].seq[consensus_start:consensus_end + 1]
          if '?' in gene_sequence:
            print ('Gene {} contains \'?\' characters indicating missing coverage. Skipping.'.format(gene))
            with open('MissingCoverage.txt', 'a') as file:
              file.write(SequenceName + '_' + gene + '\n')
            continue
          else:
            output_file.write('>' + SequenceName + '_' + gene + '\n' + str(gene_sequence) + '\n')
          
        else:
          print ('Gene {} not found in alignment.'.format(gene)) # untested output 
else:
  print ('Could not find sys.argv[2] in the provided references.') # untested output

# Check length is expected

### Take the extracted genes and reference and extract to individual files
# Is there a way to avoid doing this? But I think VIRULIGN requires a .fasta (or xml) input
# Can also just output them individually in the previous steps but having them together seemed useful
# Optionally delete temp/single gene files after? Maybe if option is switched on to keep single gene
  # files then don't name to temp, otherwise name temp_ and delete after

# Extract the sample sequences to individual fastas for VIRULIGN processing
# maybe make as temp_ files and deleted after processing
multi_fasta_sample = SequenceName + '_GenesExtracted.fasta' # python2?

for record in SeqIO.parse(multi_fasta_sample, 'fasta'):
  individual_file = 'temp_{}_only.fasta'.format(record.id)
  with open(individual_file, 'w') as file:
    file.write('>' + str(record.id) + '\n' + str(record.seq))

# Extract the reference sequences to individual fastas for VIRULIGN processing
multi_fasta_ref = SequenceName + '_Reference_GenesExtracted.fasta'

for record in SeqIO.parse(multi_fasta_ref, 'fasta'):
  individual_file = 'temp_{}_Reference_only.fasta'.format(record.id)
  with open(individual_file, 'w') as file:
    file.write('>' + str(record.id) + '\n' + str(record.seq))

"""

python2 -c "$Python_Extract_Genes" "$ReferenceAlignment" "$RefSequenceName_shell" "$GeneCoordInfo" "$GAG" "$POL" "$VIF" "$VPR" "$VPU" "$ENV" "$NEF" "$Remap"

# add warning if some are added to missing coverage

###################################################################################
# No longer needed due to mafft extraction - leaving in temporarily in case code is useful later

# Summary: 
# BLASTn full sequence to gene to find closest reference then extract gene and reference 


# ReferencePOL="/Users/s.blundell/shiver-unpaired/CodonCorrectionReferences/POL/HIV1_ALL_2021_pol_DNA.fasta"

# ### BLASTn the target sequence to a specific gene
#   # Call as function to avoid specific file names
# SampleSequence="/Users/s.blundell/codon-awareness/HIV-COM/HIV_COM_B.FR_shifted.fasta"
# blastn -query "$SampleSequence" -db ~/shiver/CodonCorrectionReferences/POL/blast_db_pol_dna -outfmt \
# "10 qseqid qseq evalue sseqid pident qlen qstart qend sstart send frames" \
# -out POL_blastn.out -max_target_seqs 1

# blast_output="POL_blastn.out"

# # Find the sseqid associated with the longest BLAST fragment (intended as the closest match)
# ClosestMatch=$(awk 'BEGIN { longest_fragment=0 } { if (length($2) > longest_fragment) { longest_fragment=length($2); \
# closest_match=$4 } } END { print closest_match }' FS=',' "$blast_output")

# echo "ClosestMatch:" "$ClosestMatch" # TESTING

# ### Extract the closest BLAST hit to a reference .fasta file
# python2 - "$ReferencePOL" "$ClosestMatch" <<END
# import sys
# from Bio import SeqIO
# from Bio.SeqRecord import SeqRecord

# ReferenceSequence = None
# with open(sys.argv[1], "r") as handle:
#   for record in SeqIO.parse(handle, "fasta"):
#     if record.id == sys.argv[2]:
#       ReferenceSequence = record.seq
#       break

# if ReferenceSequence:
#   output_file = "ReferencePOLExtracted.fasta"
#   ungapped_record = SeqRecord(id=record.id, seq=ReferenceSequence.ungap("-"), description="")
#   with open(output_file, "w") as output_handle:
#     SeqIO.write([ungapped_record], output_handle, "fasta")
# else:
#     # check how to print properly
#     print("Sequence", sys.argv[2], "was not found within the reference .fasta file.")
# END

# ### Extract the gene from the BLASTn output
# # REMOVE/skip if you already have specific genes, eg POPART data. 
# # Set up as function to avoid file paths
# # Extract qstart and qend from each line and get the min/max
# SequenceStart=$(awk -F',' '{print $7}' "$blast_output" | sort -n | head -n 1)
# SequenceEnd=$(awk -F',' '{print $8}' "$blast_output" | sort -rn | head -n 1)

# # Extract the sequence between qstart and qend, i.e. the full gene
# python2 - "$SequenceStart" "$SequenceEnd" "$SampleSequence" <<END
# import sys
# from Bio import SeqIO

# # Some key files
# sample_file = sys.argv[3]
# output_fasta = "BlastExtractedSequence.fasta"

# with open(sample_file, 'r') as blast_file:
#   first_line = blast_file.readline().strip()
#   sequence_id = first_line.split(',')[0][1:]  

# # Read the fasta file and extract the gene
# sequences_dict = SeqIO.to_dict(SeqIO.parse(sample_file, "fasta"))
# gene_sequence = sequences_dict[sequence_id][int(sys.argv[1]):int(sys.argv[2])]

# with open(output_fasta, "w") as output_handle:
#     SeqIO.write(gene_sequence, output_handle, "fasta")
# END

# SamplePOLGene="BlastExtractedSequence.fasta"

# ReferencePOLGene="ReferencePOLExtracted.fasta"
###################################################################################

### Run VIRULIGN

# Make a directory for the failed alignments
  # Add path check error messages
if [ ! -d "FailedVIRULIGN" ]; then
  mkdir FailedVIRULIGN
fi

# Determine VIRULIGN options
# Probably far better ways to structure this section - return to
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

VirulignOutput=("Nucleotides" "AminoAcids" "Mutations" "MutationTable")

if [[ "$5" == "all" ]]; then 
  Nucleotides=true
  AminoAcids=true
  Mutations=true
  MutationTable=true
else
  # Select gene output based on VirulignOPtions args
  IFS=',' read -ra listed_options <<< "$5"
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

# Print chosen options
listed_options=()
for option in "${VirulignOutput[@]}"; do
  if [ "${!option}" = true ]; then
    listed_options+=(" $option")
  fi
done

if [ "${#listed_options[@]}" -gt 0 ]; then
  echo "Outputting files for chosen options: ${listed_options[@]}"
else
  echo "Problem with virulign options. Quitting." >&2
  exit 1
fi

# Initialise variables
export_alphabet=""
export_kind=""
FailedDebug=""
FileAppend=""

# Define a function to call virulign
# Addended _func because otherwise it would presumably overwrite variables in rest of script? or can they be the same
function run_virulign {
  export_alphabet_func="$1"
  export_kind_func="$2"
  FailedDebug_func="$3"
  FileAppend_func="$4"

  for gene in "${genes[@]}"; do
    if [ "${!gene}" = true ]; then
      SampleGene=$(find . -type f -name "temp_${SequenceName_shell}_${gene}_only.fasta")
      ReferenceGene=$(find . -type f -name "temp_${RefSequenceName_shell}_${gene}_Reference_only.fasta")

      # Run VIRULIGN
      $VirulignLocation $ReferenceGene $SampleGene --exportKind $export_kind_func --exportAlphabet $export_alphabet_func --exportReferenceSequence yes --exportWithInsertions yes $FailedDebug_func > HIV_${gene}$FileAppend_func
    
    fi
  done
}

### Define virulign command options based on chosen output options
if [[ "$Nucleotides" == "true" ]]; then
  export_alphabet="Nucleotides"
  export_kind="GlobalAlignment"
  FailedDebug="--nt-debug ${gene}FailedVIRULIGN"
  # Currently only using this for Nucleotides. Could make diff file names for each failed output but need to consider how useful this is
  FileAppend="_Nucl_corrected.fasta"

  run_virulign "$export_alphabet" "$export_kind" "$FailedDebug" "$FileAppend" || { echo "Problem running virulign. Quitting." >&2; exit 1; }
fi

if [[ "$AminoAcids" == "true" ]]; then
  export_alphabet="AminoAcids"
  export_kind="GlobalAlignment"
  FailedDebug=""
  FileAppend="_AminoAcids_corrected.fasta"

  run_virulign "$export_alphabet" "$export_kind" "$FailedDebug" "$FileAppend" || { echo "Problem running virulign. Quitting." >&2; exit 1; }
fi

if [[ "$Mutations" == "true" ]]; then
  export_alphabet="Nucleotides"
  export_kind="Mutations"
  FailedDebug=""
  FileAppend="_Mutations.csv"

  run_virulign "$export_alphabet" "$export_kind" "$FailedDebug" "$FileAppend" || { echo "Problem running virulign. Quitting." >&2; exit 1; }
fi

if [[ "$MutationTable" == "true" ]]; then
  export_alphabet="AminoAcids"
  export_kind="MutationTable"
  FailedDebug=""
  FileAppend="_MutationTable.csv"

  run_virulign "$export_alphabet" "$export_kind" "$FailedDebug" "$FileAppend" || { echo "Problem running virulign. Quitting." >&2; exit 1; }
fi

# Check options are set correctly
if [[ -z $export_alphabet || -z $export_kind || -z $FileAppend ]]; then
    echo "Error: One or more options are not set. Quitting."
    exit 1
fi


# Print frameshift info to file
# can only happen if mutations is enabled
# TODO if sequence_name already in file return error
# Could be set up to return more meaningful info
if [[ ! -f "Frameshifts.txt"  ]] && [[ "$Mutations" == "true" ]]; then
  echo "Sequence Name, Frameshift number, Mutations" > Frameshifts.txt
fi

if [[ -f "Frameshifts.txt"  ]] && [[ "$Mutations" == "true" ]]; then
  # Extract frameshift information
  for gene in "${genes[@]}"; do
    if [[ -f "HIV_${gene}_Mutations.csv" ]]; then
      frameshift_info=$(awk -F ',' 'NR==3 {print $1 ", " $4 ", " $7}' "HIV_${gene}_Mutations.csv")
      frameshift_number=$(echo "$frameshift_info" | awk -F', ' '{print $2}')
      if [[ "$frameshift_number" != "0" ]]; then
        echo "$frameshift_info" >> Frameshifts.txt
      fi
    fi
  done
elif [[ "$Mutations" == "false" ]]; then
  echo "To list frameshifts found in the sample enable the 'Mutations' virulign option."
fi

# Delete Failed directory if empty
# Need to test with one that I know fails - maybe manually change reference. 
if [[ -d FailedVIRULIGN ]]; then
  rmdir ./FailedVIRULIGN || echo "Some sequences failed to align correctly, check /FailedVIRULIGN."
fi


# For virulign:
# Current setup adds insertions to the reference where new codons are added in the sample
# nt-debug is only added to nucleotide sequence atm, could set it up so it's done for at least one
  # of the commands, eg if only doing AA or mutations. Can have different folder names to do multiple

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

# set up help

cd ../ # Return to ~/shiver. temp for testing
