#!/usr/bin/env python

import sys
import os
from Bio import AlignIO
from Bio import SeqIO
import argparse

HelpMessage = "Script containing all the required python functions for Codon_Correction.sh (TODO)" # expand this

# Set up arguments
parser = argparse.ArgumentParser(description=HelpMessage)
parser.add_argument('FunctionName')
parser.add_argument('--X') # example optional arg - delete later. each fucntion will have it's own optional args
parser.add_argument('--InitDir', help='The location chosen for the init directory')
parser.add_argument('--GeneCoordInfo', help='Gene coordinate information file') # be more descriptive of format
parser.add_argument('--GenomeFile', help='Fasta file of all the references as complete genomes (as opposed to individual genes).')
parser.add_argument('--AlignmentFile', help='Fasta file containing the HXB2 sequence included with shiver aligned to the sample.')
parser.add_argument('--SingleSeq', help='Set to true or fasle in CC.sh depending on number of sequences in the sample file.')
parser.add_argument('--ReferenceName', help='Name of the reference that BLASTs to the gene')
parser.add_argument('--Gene', help='The current gene being processed')
parser.add_argument('--AlignmentToRef', help='Fasta file containing the reference sequence for each gene aligned to the sample genome.')
# can add defaults if useful

args = parser.parse_args()
# Select function to call
function_name = args.FunctionName

# MakeReferenceDatabase uses the whole genome reference file to  extract individual genes to a separate fasta file based on gene type
def MakeReferenceDatabase(InitDir, GeneCoordInfo, GenomeFile):
  print ('Now extracting reference genes from database')
  # Parse the whole genome sequences
  with open(GenomeFile, 'r') as genome_handle:
    genome_records = SeqIO.to_dict(SeqIO.parse(genome_handle, "fasta"))

  # Assign gene info per sequence in coordinates file
  gene_loci = {}
  with open(GeneCoordInfo, 'r') as file:
    next(file) # Skip first line containing headers
    for line in file:
      gene_info = line.strip().split(',')
      CoordRefName = gene_info[0]      
      gene_loci[CoordRefName] = {
        # Gene: (gene_start, gene_end)
        'GAG': (int(gene_info[1]), int(gene_info[2])),
        'POL': (int(gene_info[3]), int(gene_info[4])),
        'VIF': (int(gene_info[5]), int(gene_info[6])),
        'VPR': (int(gene_info[7]), int(gene_info[8])),
        'VPU': (int(gene_info[9]), int(gene_info[10])),
        'ENV': (int(gene_info[11]), int(gene_info[12])),
        'NEF': (int(gene_info[13]), int(gene_info[14]))
      }
      # Find matching genome sequence
      if CoordRefName in genome_records:
        GenomeSeq = genome_records[CoordRefName].seq
      else:
        print ('Sequence {} not found in genome file. Quitting'.format(CoordRefName))
        return 1

      # Extract reference sequence for each gene
      for gene_class, (gene_start, gene_end) in gene_loci[CoordRefName].items():
        OutputFile = os.path.join(InitDir, 'ReferenceGenes_{}.fasta'.format(gene_class))
        with open(OutputFile, 'a') as output_file:
          if gene_start is not None and gene_end is not None:
            ref_pos = 0
            ref_start = 0
            ref_end = 0
            for i, base in enumerate(GenomeSeq):
              if base != '-':
                ref_pos += 1
              if ref_pos == gene_start:
                ref_start = i
              if ref_pos == gene_end:
                ref_end = i
                gene_sequence = GenomeSeq[ref_start:ref_end + 1].ungap('-')
                output_file.write('>' + CoordRefName + '_' + gene_class + '\n' + str(gene_sequence) + '\n')    
                break
  return 0

# Call MakeReferenceDatabase
if function_name == 'MakeReferenceDatabase':
  MakeReferenceDatabase(args.InitDir, args.GeneCoordInfo, args.GenomeFile)

# ExtractWithHXB2 uses a set of predetermined HXB2 coordinates to extract genes from the sample file. Using a set reference
# will reduce accuracy of results, but the extracted genes are only used to BLASTn to a closer reference which will then
# be used to extract the final sample gene sequences used for VIURLIGN
def ExtractWithHXB2(AlignmentFile, SingleSeq):

  SingleSequence_check = SingleSeq
  RefSequenceNumber = 0
  if SingleSequence_check == 'true':
    RefSequenceNumber = 1
  else:
    RefSequenceNumber = 2

  # Assign variables
  HXB2Alignment = AlignIO.read(AlignmentFile, 'fasta')
  SequenceName = HXB2Alignment[0].id
  RefSequenceLength = len(HXB2Alignment[RefSequenceNumber])
  output_file_name = 'temp_preBLAST_Sample_Genes.fasta' # NEEDS TO BE PER GENE

  # Extract gene coordinates from the reference file - for HXB2 these are set to predetermined values
  hxb2_loci = {
    # Gene: (gene_start, gene_end)
    'GAG': (int(790), int(2292)),
    'POL': (int(2253), int(5096)),
    'VIF': (int(5041), int(5619)),
    'VPR': (int(5559), int(5850)),
    'VPU': (int(6062), int(6310)),
    'ENV': (int(6225), int(8795)),
    'NEF': (int(8797), int(9168))
  }

  # Gene extraction for sample sequence
  if hxb2_loci:
    for gene, (gene_start, gene_end) in hxb2_loci.items():
      OutputFile = "temp_preBLAST_Sample_Gene_{}.fasta".format(gene) # is this correct?
      with open(OutputFile, 'w') as output_file:
        if gene_start is not None and gene_end is not None:
          # Iterate over the sequence to determine gene coordinates of the sample sequence
          consensus_start = 0
          consensus_end = 0
          sequence_pos = 0
          for i in xrange(RefSequenceLength): # try ennumerate as above, and/or set sequence to a variable
            if HXB2Alignment[RefSequenceNumber][i] != '-':
              sequence_pos += 1
            if sequence_pos == gene_start:
              consensus_start = i
            if sequence_pos == gene_end: 
              consensus_end = i
              break

          # Write the gene sequence to file unless it contains '?'
          if consensus_start is not None and consensus_end is not None:
            gene_sequence = HXB2Alignment[0].seq[consensus_start:consensus_end + 1].ungap('-')
            if '?' in gene_sequence:
              print ('Gene {} contains \'?\' characters indicating missing coverage. Skipping.'.format(gene))
              with open('MissingCoverage.txt', 'a') as file:
                file.write(SequenceName + '_' + gene + '\n')
              continue
            else:
              output_file.write('>' + SequenceName + '_' + gene + '\n' + str(gene_sequence) + '\n')
  return 0


# Call ExtractWithHXB2
if function_name == 'ExtractWithHXB2':
  ExtractWithHXB2(args.AlignmentFile, args.SingleSeq)

# ExtractRefFromFasta takes the single whole genome sequence out of the full reference fasta
def ExtractRefFromFasta(ReferenceName, GenomeFile, Gene):
  output_ref = 'temp_Reference_Genome_' + Gene + '.fasta' 
  found = False
  with open(output_ref, 'w') as file:
    for record in SeqIO.parse(GenomeFile, 'fasta'):
      if record.id == ReferenceName:
        SeqIO.write(record, file, 'fasta')
        found = True
        break

  if not found:
    print('Reference sequence {} was not found win the provided reference file {}.'.format(ReferenceName, GenomeFile))
  return 0

# Call ExtractRefFromFasta
if function_name == 'ExtractRefFromFasta':
  ExtractRefFromFasta(args.ReferenceName, args.GenomeFile, args.Gene)

def ExtractRefandSample(AlignmentToRef, SingleSeq, ReferenceName, GeneCoordInfo, Gene):

  SingleSequence_check = SingleSeq
  RefSequenceNumber = 0
  if SingleSequence_check == 'true':
    RefSequenceNumber = 1
  else:
    RefSequenceNumber = 2

  # Assign variables
  ReferenceAlignment = AlignIO.read(AlignmentToRef, 'fasta')
  SequenceName = ReferenceAlignment[0].id
  RefSequenceLength = len(ReferenceAlignment[RefSequenceNumber])
  output_file_name = '{}_GenesExtracted.fasta'.format(SequenceName)
  ref_output_file_name = '{}_Reference_GenesExtracted.fasta'.format(SequenceName) # Could be changed to not use SequenceName if that's clearer

  # Extract gene coordinates from the reference file
  with open(GeneCoordInfo, 'r') as file:
    for line in file:
      gene_info = line.strip().split(',')
      if gene_info[0] == ReferenceName:
        # Identify gene and select corresponding gene info coordinates
        if Gene == "GAG":
          gene_start = int(gene_info[1])
          gene_end = int(gene_info[2])
        elif Gene == "POL":
          gene_start = int(gene_info[3])
          gene_end = int(gene_info[4])
        elif Gene == "VIF":
          gene_start = int(gene_info[5])
          gene_end = int(gene_info[6])
        elif Gene == "VPR":
          gene_start = int(gene_info[7])
          gene_end = int(gene_info[8])
        elif Gene == "VPU":
          gene_start = int(gene_info[9])
          gene_end = int(gene_info[10])
        elif Gene == "ENV":
          gene_start = int(gene_info[11])
          gene_end = int(gene_info[12])
        elif Gene == "NEF":
          gene_start = int(gene_info[13])
          gene_end = int(gene_info[14])
        else:
          print ('Could not match gene name to any of the expected inputs. Quitting')
          return 1
        break
      else:
        continue
    if gene_start is None or gene_end is None:
      print ('Could not match the reference sequence to any within the reference file. Quitting.')
      return 1

  # Gene extraction for reference sequence
  with open(ref_output_file_name, 'a') as output_file:
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
      # Write the reference sequence to a file
      reference_sequence = ReferenceAlignment[RefSequenceNumber].seq[ref_start:ref_end + 1].ungap('-')
      output_file.write('>' + ReferenceAlignment[RefSequenceNumber].id + '_' + Gene + '\n' + str(reference_sequence) + '\n')
    else:
      print ('Unable to extract the reference coordinates. Quitting.')
      return 1

    # Check length of ref is as expected TODO

  # Gene extraction for sample sequence
  with open(output_file_name, 'a') as output_file:
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

    if consensus_start is not None and consensus_end is not None:
      # Calculate indel size
      gene_sequence_gapped = ReferenceAlignment[0].seq[consensus_start:consensus_end + 1]
      sample_gaps = gene_sequence_gapped.count('-')
      indel_difference =  (consensus_end + 1 - consensus_start - sample_gaps) - (gene_end + 1 - gene_start)
      if indel_difference != 0:
        print ('\033[93mChange in length of {} in {}_{} relative to the reference\033[0m'.format(indel_difference, SequenceName, Gene))

      # Write the gene sequence to file unless it contains '?'
      gene_sequence = ReferenceAlignment[0].seq[consensus_start:consensus_end + 1].ungap('-')
      if '?' in gene_sequence:
        print ('Gene {} contains \'?\' characters indicating missing coverage. Skipping.'.format(gene))
        with open('MissingCoverage.txt', 'a') as file:
          file.write(SequenceName + '_' + Gene + '\n')
      else:
        output_file.write('>' + SequenceName + '_' + Gene + '\n' + str(gene_sequence) + '\n')
      
    else:
      print ('Unable to extract sample coordinates for {}.'.format(Gene)) # untested output 

  # Check length is expected

  # Extract the sample sequences to individual fastas for VIRULIGN processing
  multi_fasta_sample = SequenceName + '_GenesExtracted.fasta'

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

  # Return error if the [Sample]_GenesExtracted.fasta file is empty
  if os.path.getsize(multi_fasta_sample) == 0:
    print('{} is empty. This can happen if all genes being analysed have missing coverage (? characters). Check MissingCoverage.txt. Now quitting.'.format(multi_fasta_sample))
    return 1


# Call ExtractRefandSample
if function_name == 'ExtractRefandSample':
  ExtractRefandSample(args.AlignmentToRef, args.SingleSeq, args.ReferenceName, args.GeneCoordInfo, args.Gene)

