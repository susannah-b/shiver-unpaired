#!/usr/bin/env python

import sys
import os
from Bio import AlignIO
from Bio import SeqIO

# This script functions within CodonCorrection.sh, to extract the specified genes from both the provided sample file
# and reference. 

# Determine position of reference sequence in alignment file
# add check for RefSequenceNumber is not 0 later

# add any necessary checks before running

SingleSequence_check = sys.argv[11]
RefSequenceNumber = 0
if SingleSequence_check == 'true':
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

# Gene extraction for reference sequence
gene_reference_gaps = {} # part of my indel detection for testing - remove reference gaps calcs for final code as it's not foolproof (if gaps are in the last section equal to the # of gaps)
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
            # Count gaps for indel calculation
            reference_sequence_gapped = ReferenceAlignment[RefSequenceNumber].seq[ref_start:ref_end + 1] # if this part ignored gaps for ref end determination it would work for all sequences
            gene_reference_gaps[gene] = reference_sequence_gapped.count('-')

            # Write the reference sequence to a file
            reference_sequence = ReferenceAlignment[RefSequenceNumber].seq[ref_start:ref_end + 1].ungap('-')
            output_file.write('>' + ReferenceAlignment[RefSequenceNumber].id + '_' + gene + '\n' + str(reference_sequence) + '\n')

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

        if consensus_start is not None and consensus_end is not None:
          # Calculate indel size
          gene_sequence_gapped = ReferenceAlignment[0].seq[consensus_start:consensus_end + 1]
          sample_gaps = gene_sequence_gapped.count('-')
          indel_difference =  (consensus_end + 1 - consensus_start - sample_gaps) - (gene_end + 1 - gene_start)
          if indel_difference != 0:
            print ('Change in length of {} in {}_{} relative to the reference'.format(indel_difference, SequenceName, gene))

          # Write the gene sequence to file unless it contains '?'
          gene_sequence = ReferenceAlignment[0].seq[consensus_start:consensus_end + 1].ungap('-')
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
  sys.exit(1)