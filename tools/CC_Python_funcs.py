#!/usr/bin/env python

import sys
import os
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
# can add defaults if useful

args = parser.parse_args()
# Select function to call
function_name = args.FunctionName

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






