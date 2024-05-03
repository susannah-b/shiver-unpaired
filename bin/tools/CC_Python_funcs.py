#!/usr/bin/env python

import sys
import os
from Bio import AlignIO
from Bio import SeqIO
import csv
import argparse
from AuxiliaryFunctions import ungap

HelpMessage = "Script containing all the required python functions for Codon_Correction.sh. Each function is called with the relevant \
optional arguments in the CodonCorrection script and init script."

# Set up arguments
parser = argparse.ArgumentParser(description=HelpMessage)
parser.add_argument('FunctionName')
parser.add_argument('--InitDir', help='The location chosen for the init directory')
parser.add_argument('--GeneCoordInfo', help='Gene coordinate information file') # be more descriptive of format
parser.add_argument('--GenomeFile', help='Fasta file of all the references as complete genomes (as opposed to individual genes).')
parser.add_argument('--AlignmentFile', help='Fasta file containing the HXB2 sequence included with shiver aligned to the sample.')
parser.add_argument('--SingleSeq', help='Set to true or false in CC.sh depending on number of sequences in the sample file.')
parser.add_argument('--ReferenceName', help='Name of the reference that BLASTs to the gene')
parser.add_argument('--Gene', help='The current gene being processed')
parser.add_argument('--AlignmentToRef', help='Fasta file containing the reference sequence for each gene aligned to the sample genome.')
parser.add_argument('--OutputFile', help='The name of the output fasta file to write to')
parser.add_argument('--InputFile', help='The name of the input fasta file')
parser.add_argument('--SequenceNumber', help='The number of the sequence needed to extract (1-based)')
parser.add_argument('--OutputCSV', help='The name of the output .CSV file to write to')
parser.add_argument('--SequenceName', help='Name of the sample sequence')
parser.add_argument('--ExtractedGeneFile', help='The name of the sample gene file extracted after alignment with the BLAST reference')
args = parser.parse_args()

# Select function to call
function_name = args.FunctionName

# MakeReferenceDatabase uses the whole genome reference file to  extract individual genes to a separate fasta file based on gene type
def MakeReferenceDatabase(InitDir, GeneCoordInfo, GenomeFile):
  print('Now extracting reference genes from database')
  # Parse the whole genome sequences
  try:
    with open(GenomeFile, 'r') as genome_handle:
      genome_records = SeqIO.to_dict(SeqIO.parse(genome_handle, "fasta"))
  except IOError:
    print(('Error: Could not open file {}'.format(GenomeFile)))
    raise
  except ValueError:
    print(('Error: File {} is not in the expected format'.format(GenomeFile)))
    raise

  # Assign gene info per sequence in coordinates file
  gene_loci = {}
  try:
    with open(GeneCoordInfo, 'r') as file:
      next(file)
      for line in file:
        gene_info = line.strip().split(',')
        if len(gene_info) != 15:
          print('Error: The gene coordinate file has an unexpected number of fields. Quitting.')
          exit(1)
        
        CoordRefName = gene_info[0]
        if CoordRefName in gene_loci:
          print('Error: encountered reference', CoordRefName,
          'a second time in', GeneCoordInfo, file=sys.stderr)
          exit(1)
        
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
        
        # Check start & end coords are OK
        for gene, (start, end) in gene_loci[CoordRefName].items():
          if end < start:
            print("Error: for reference", CoordRefName + ", gene", gene,
            "has an end position before its start position.", file=sys.stderr)
            exit(1)
          if (end - start + 1) % 3 != 0:
            print("Error: for reference", CoordRefName + ", gene", gene,
            "has length", end - start + 1, "which is not a multiple of 3. This",
            "is not compatible with use for codon correction.",
            file=sys.stderr)
            exit(1)
        
        # Find matching genome sequence
        if CoordRefName in genome_records:
          GenomeSeq = genome_records[CoordRefName].seq
        else:
          print(('Sequence {} not found in genome file. Quitting'.format(CoordRefName)))
          exit(1)

        # Extract reference sequence for each gene
        for gene_class, (gene_start, gene_end) in list(gene_loci[CoordRefName].items()):
          OutputFile = os.path.join(InitDir, 'ReferenceGenes_{}.fasta'.format(gene_class))
          try:
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
                    gene_sequence = ungap(GenomeSeq[ref_start:ref_end + 1])
                    output_file.write('>' + CoordRefName + '_' + gene_class + '\n' + str(gene_sequence) + '\n')    
                    break
          except IOError:
            print(('Error: Could not open file {}'.format(OutputFile)))
            raise
          except ValueError:
            print(('Error: File {} is not in the expected format'.format(OutputFile)))
            raise
    missing_refs = [ref_ for ref_ in genome_records.keys() if not ref_ in \
    gene_loci.keys()]
    if len(missing_refs) > 0:
      print("Warning: some references have a sequence in ", GenomeFile,
      " but their gene coordinates are not given in ", GeneCoordInfo, 
      ". They will not be used for codon correction. These references are: ",
      ' '.join(missing_refs), sep='', file=sys.stderr)
  except IOError:
    print(('Error: Could not open file {}'.format(GeneCoordInfo)))
    raise
  except ValueError:
    print(('Error: File {} is not in the expected format'.format(GeneCoordInfo)))
    raise
  return 0

# Call MakeReferenceDatabase
if function_name == 'MakeReferenceDatabase':
  MakeReferenceDatabase(args.InitDir, args.GeneCoordInfo, args.GenomeFile)

# ExtractWithHXB2 uses a set of predetermined HXB2 coordinates to extract genes from the sample file. Using a set reference
# will reduce accuracy of results, but the extracted genes are only used to BLASTn to a closer reference which will then
# be used to extract the final sample gene sequences used for VIRULIGN
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
    for gene, (gene_start, gene_end) in list(hxb2_loci.items()):
      OutputFile = "temp_preBLAST_Sample_Gene_{}.fasta".format(gene) # is this correct?
      try:
        with open(OutputFile, 'w') as output_file:
          if gene_start is not None and gene_end is not None:
            # Iterate over the sequence to determine gene coordinates of the sample sequence
            consensus_start = 0
            consensus_end = 0
            sequence_pos = 0
            for i in range(RefSequenceLength): # try ennumerate as above, and/or set sequence to a variable
              if HXB2Alignment[RefSequenceNumber][i] != '-':
                sequence_pos += 1
              if sequence_pos == gene_start:
                consensus_start = i
              if sequence_pos == gene_end: 
                consensus_end = i
                break

            # Write the gene sequence to file unless it contains '?'
            if consensus_start is not None and consensus_end is not None:
              gene_sequence = ungap(HXB2Alignment[0].seq[consensus_start:consensus_end + 1])
              if '?' in gene_sequence:
                print(('Gene {} contains \'?\' characters indicating missing coverage. Skipping.'.format(gene)))
                with open('MissingCoverage.txt', 'a') as file:
                  file.write(SequenceName + '_' + gene + '\n')
                continue
              else:
                output_file.write('>' + SequenceName + '_' + gene + '\n' + str(gene_sequence) + '\n')
      except IOError:
        print(('Error: Could not open file {}'.format(OutputFile)))
        raise
      except ValueError:
        print(('Error: File {} is not in the expected format'.format(OutputFile)))
        raise
        
  return 0


# Call ExtractWithHXB2
if function_name == 'ExtractWithHXB2':
  ExtractWithHXB2(args.AlignmentFile, args.SingleSeq)

# ExtractRefFromFasta takes the single whole genome sequence out of the full reference fasta
def ExtractRefFromFasta(ReferenceName, GenomeFile, Gene):
  output_ref = 'temp_Reference_Genome_' + Gene + '.fasta' 
  found = False
  try:
    with open(output_ref, 'w') as file:
      for record in SeqIO.parse(GenomeFile, 'fasta'):
        if record.id == ReferenceName:
          SeqIO.write(record, file, 'fasta')
          found = True
          break
  except IOError:
    print(('Error: Could not open file {}'.format(GenomeFile)))
    raise
  except ValueError:
    print(('Error: File {} is not in the expected format'.format(GenomeFile)))
    raise

  if not found:
    print(('Reference sequence {} was not found within the provided reference file {}.'.format(ReferenceName, GenomeFile)))
  
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
  ref_output_file_name = '{}_Reference_GenesExtracted.fasta'.format(SequenceName) # Could be changed to use ReferenceName if that's clearer

  # Extract gene coordinates from the reference file
  try:
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
            print('Could not match gene name to any of the expected inputs. Quitting')
            exit(1)
          break
        else:
          continue
      if gene_start is None or gene_end is None:
        print('Could not match the reference sequence to any within the reference file. Quitting.')
        exit(1)
  except IOError:
    print(('Error: Could not open file {}'.format(GeneCoordInfo)))
    raise
  except ValueError:
    print(('Error: File {} is not in the expected format'.format(GeneCoordInfo)))
    raise

  # Gene extraction for reference sequence
  try:
    with open(ref_output_file_name, 'a') as output_file:
      ref_pos = 0
      ref_start = 0
      ref_end = 0
      for i in range(RefSequenceLength):
        if ReferenceAlignment[RefSequenceNumber][i] != '-':
          ref_pos += 1
        if ref_pos == gene_start:
          ref_start = i
        if ref_pos == gene_end:
          ref_end = i
          break

      if ref_start is not None and ref_end is not None:
        # Write the reference sequence to a file
        reference_sequence = ungap(ReferenceAlignment[RefSequenceNumber].seq[ref_start:ref_end + 1])
        output_file.write('>' + ReferenceAlignment[RefSequenceNumber].id + '_' + Gene + '\n' + str(reference_sequence) + '\n')
      else:
        print('Unable to extract the reference coordinates. Quitting.')
        exit(1)
  except IOError:
    print(('Error: Could not open file {}'.format(ref_output_file_name)))
    raise
  except ValueError:
    print(('Error: File {} is not in the expected format'.format(ref_output_file_name)))
    raise

  # Gene extraction for sample sequence
  try:
    with open(output_file_name, 'a') as output_file:
      # Iterate over the sequence to determine gene coordinates of the sample sequence
      consensus_start = None
      consensus_end = None
      sequence_pos = 0
      indel_difference = 0
      for i in range(RefSequenceLength):
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
          print(('Change in length of {} in {}_{} relative to the reference'.format(indel_difference, SequenceName, Gene))) # This can be removed if not useful, it was mostly used when I was debugging length issues

        # Write the gene sequence to file unless it contains '?'
        gene_sequence = ungap(ReferenceAlignment[0].seq[consensus_start:consensus_end + 1])
        if '?' in gene_sequence:
          print(('Gene {} contains \'?\' characters indicating missing coverage. Skipping.'.format(Gene)))
          with open('MissingCoverage.txt', 'a') as file:
            file.write(SequenceName + '_' + Gene + '\n')
        else:
          output_file.write('>' + SequenceName + '_' + Gene + '\n' + str(gene_sequence) + '\n')

        # Write the sample gene start to a file for later use in indel position calculation
        gene_start_full_genome = len(ungap(ReferenceAlignment[0].seq[0:consensus_start + 1]))
        with open('temp_GeneStarts.csv', 'a', newline='') as csvfile:
          fieldnames = ['Sequence', 'Gene', 'GeneStart']
          writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
          if csvfile.tell() == 0:
            writer.writeheader()
          writer.writerow({'Sequence': SequenceName, 'Gene': Gene, 'GeneStart': gene_start_full_genome})
      else:
        print(('Unable to extract sample coordinates for {}.'.format(Gene))) 
  except IOError:
    print(('Error: Could not open file {}'.format(output_file_name)))
    raise
  except ValueError:
    print(('Error: File {} is not in the expected format'.format(output_file_name)))
    raise

  # Extract the sample sequences to individual fastas for VIRULIGN processing
  multi_fasta_sample = SequenceName + '_GenesExtracted.fasta'

  for record in SeqIO.parse(multi_fasta_sample, 'fasta'):
    individual_file = 'temp_{}_only.fasta'.format(record.id)
    try:
      with open(individual_file, 'w') as file:
        file.write('>' + str(record.id) + '\n' + str(record.seq))
    except IOError:
      print(('Error: Could not open file {}'.format(individual_file)))
      raise
    except ValueError:
      print(('Error: File {} is not in the expected format'.format(individual_file)))
      raise

  # Extract the reference sequences to individual fastas for VIRULIGN processing
  multi_fasta_ref = SequenceName + '_Reference_GenesExtracted.fasta'

  for record in SeqIO.parse(multi_fasta_ref, 'fasta'):
    individual_file = 'temp_{}_Reference_only.fasta'.format(record.id)
    try:
      with open(individual_file, 'w') as file:
        file.write('>' + str(record.id) + '\n' + str(record.seq))
    except IOError:
      print(('Error: Could not open file {}'.format(individual_file)))
      raise
    except ValueError:
      print(('Error: File {} is not in the expected format'.format(individual_file)))

  # Return error if the [Sample]_GenesExtracted.fasta file is empty
  if os.path.getsize(multi_fasta_sample) == 0:
    print(('{} is empty. This can happen if all genes being analysed have missing coverage (? characters). Check MissingCoverage.txt. Now quitting.'.format(multi_fasta_sample)))
    exit(1)
  
  return 0

# Call ExtractRefandSample
if function_name == 'ExtractRefandSample':
  ExtractRefandSample(args.AlignmentToRef, args.SingleSeq, args.ReferenceName, args.GeneCoordInfo, args.Gene)

# ExtractSequence takes the specified sequence out of the specified file
def ExtractSequence(OutputFile, InputFile, SequenceNumber):
  try:
    with open(OutputFile, 'w') as file:
      InputFasta=AlignIO.read(InputFile, 'fasta')
      Sequence=InputFasta[int(SequenceNumber) - 1]
      SeqIO.write(Sequence, file, 'fasta')
  except IOError:
    print(('Error: Could not open file {}'.format(InputFile)))
    raise
  except ValueError:
    print(('Error: File {} is not in the expected format'.format(InputFile)))
    raise
  return 0

# Call ExtractSequence
if function_name == 'ExtractSequence':
  ExtractSequence(args.OutputFile, args.InputFile, args.SequenceNumber)


# CategoriseIndels classfiies indels as either insertions or deletions based on gaps in the reference
#  This is slightly biased towards insertions as it will classify an indel as insertion in corrected regions that are not present in the reference (hence gapped)
def CategoriseIndels(InputFile, OutputCSV, Gene, SequenceName, ReferenceName, GeneCoordInfo, GenomeFile, AlignmentFile, ExtractedGeneFile):
  SampleNs = []
  indels = {}
  try:
    # Find pre-existing N bases in the sample gene
    with open(ExtractedGeneFile, 'r') as file:
      ExtractedGene = SeqIO.read(file, 'fasta')
      ExtractedSeq = ExtractedGene.seq
      if 'N' in ExtractedSeq or 'n' in ExtractedSeq:
        sequence_pos = 0
        for i, base in enumerate(ExtractedSeq):
          if base != '-':
            sequence_pos += 1
          if base == 'N' or base == 'n':
            SampleNs.append(int(sequence_pos))
  except IOError:
    print(('Error: Could not open file {}'.format(ExtractedGeneFile)))
    raise
  except ValueError:
    print(('Error: File {} is not in the expected format'.format(ExtractedGeneFile)))
    raise
      
  try:
    # Determine indel type and position
    with open(InputFile, 'r') as file:
      sequences = list(SeqIO.parse(file, 'fasta'))
      reference = sequences[0].seq
      corrected = sequences[1].seq
      IndelCount = 0
      NCount = 0
      IndelModifier = 0
      sample_pos = 0
      gap_count = 0
      for i, base in enumerate(corrected):
        if base != '-':
          sample_pos += 1
        else:
          gap_count += 1
        if base == 'N':
          if sample_pos in SampleNs: # TODO set this up so if sampleNs is empty it won't compare (ie most sequences)
            continue
          else:
            # Count added N bases in order to determine indel position relative to the original gene
            NCount += 1
            # Increase downstream SampleNs by 1 to account for added bases relative to the original gene extraction
            SampleNs = [x + 1 if x > sample_pos else x for x in SampleNs]
            # Check if the prior base is not a VIRULIGN-added N - in which case it is included in the prior indel (so only the first base of each indel is used)
            priorbase = sample_pos - 1
            if corrected[i - 1] != 'N' or priorbase in SampleNs:
              if reference[i] == '-':
                # Modifier to give the position of the actual indel as opposed to the position of the virulign-added Ns. Uses position in the codon
                # to determine how virulign has corrected the frameshift - i.e. if the N correction is before or after the frameshifting base
                if (i - gap_count) % 3 == 0:
                  IndelModifier = 0 # (first base of the codon is N so indel pos needs no adjustment - eg NNB/NBB)
                elif (i - gap_count) % 3 == 1:
                  IndelModifier = 1 # (position 1 is an inserted base so indel pos needs to be adjusted by 1 - eg BNN/BNB)
                else:
                  if (i - 2) == 'N' and (sample_pos -2) not in SampleNs:
                    continue # (if NBN then the third N should not be classified as a separate indel)
                  else:
                    IndelModifier = 2 # (first two bases are inserted so indel pos needs to be adjusted by 2 - i.e. BBN)
                # Position of current indel with previous Ns subtracted and -1 to record insertion eith respect to last 'real' base
                indel_pos = int(sample_pos) - (int(NCount) - 1) - 1 - IndelModifier
                indels["Indel{}".format(i)] = {"Type": "Insertion", "GenePosition": indel_pos, "Gene": Gene, "Sequence": SequenceName} # Note that IndelN is i value not base position
                IndelCount += 1
              else:
                # Position of current indel with previous Ns subtracted and -1 to record deletion with respect to last 'real' base
                indel_pos = int(sample_pos) - (int(NCount) -1) - 1
                indels["Indel{}".format(i)] = {"Type": "Deletion", "GenePosition": indel_pos, "Gene": Gene, "Sequence": SequenceName}
                IndelCount += 1
              # Determine indel position relative to the full genome
              try: 
                with open('temp_GeneStarts.csv', 'r') as csvfile:
                  reader = csv.DictReader(csvfile)
                  for row in reader:
                    if row['Sequence'] == SequenceName and row['Gene'] == Gene:
                      gene_start_pos = int(row['GeneStart'])
                      genome_indel_pos = int(gene_start_pos) - 1 + int(indel_pos)
                  indels["Indel{}".format(i)].update({"GenomePosition": genome_indel_pos})
              except IOError:
                print(('Error: Could not open file {}'.format('temp_GeneStarts.csv')))
                raise
              except ValueError:
                print(('Error: File {} is not in the expected format'.format('temp_GeneStarts.csv')))
                raise  
              # Determine indel position relative to the full reference sample alignment - based on the last 'real' ref base so where gaps occur in the ref the position may be imprecise
              CorrectedAlignment = AlignIO.read(InputFile, 'fasta')
              ref_indel_pos = len(ungap(CorrectedAlignment[0].seq[0:i])) # For the gene only, not genome
              ref_gene_start_pos = None
              try:
                with open(GeneCoordInfo, 'r') as file:
                  for line in file:
                    gene_info = line.strip().split(',')
                    if gene_info[0] == ReferenceName:
                      # Assign reference gene starts
                      if Gene == "GAG":
                        ref_gene_start_pos = gene_info[1]
                      elif Gene == "POL":
                        ref_gene_start_pos = gene_info[3]
                      elif Gene == "VIF":
                        ref_gene_start_pos = gene_info[5]
                      elif Gene == "VPR":
                        ref_gene_start_pos = gene_info[7]
                      elif Gene == "VPU":
                        ref_gene_start_pos = gene_info[9]
                      elif Gene == "ENV":
                        ref_gene_start_pos = gene_info[11]
                      elif Gene == "NEF":
                        ref_gene_start_pos = gene_info[13]
                      else:
                        print('Could not match gene name to any of the expected inputs. Quitting')
                        exit(1)
                      break
                    else:
                      continue
                  if ref_gene_start_pos is None:
                    print('Could not match the reference sequence to any within the reference file. Quitting.')
                    exit(1)
              except IOError:
                print(('Error: Could not open file {}'.format(GeneCoordInfo)))
                raise
              except ValueError:
                print(('Error: File {} is not in the expected format'.format(GeneCoordInfo)))
                raise
              if ref_gene_start_pos is not None:
                ref_genome_indel_pos = int(ref_gene_start_pos) - 1 + int(ref_indel_pos) # Not included in Frameshifts.csv
                try: 
                  with open(GenomeFile, 'r') as file:
                    genome_records = SeqIO.to_dict(SeqIO.parse(file, "fasta"))
                    for record_id, record in list(genome_records.items()):
                      if record_id == ReferenceName:
                        GlobalSeq = record.seq
                        sequence_pos = 0
                        for j in range(len(GlobalSeq)):
                          if GlobalSeq[j] != '-':
                            sequence_pos += 1
                          if sequence_pos == ref_genome_indel_pos:
                            global_ref_indel_pos = int(j) + 1
                            indels["Indel{}".format(i)].update({"GlobalGenomeRefPosition": global_ref_indel_pos})
                            break

                except IOError:
                  print(('Error: Could not open file {}'.format(GenomeFile)))
                  raise
                except ValueError:
                  print(('Error: File {} is not in the expected format'.format(GenomeFile)))
                  raise
      
              else:
                print(('Unable to find the reference gene start for {}. Quitting.'.format(Gene)))
              
              # Determine position relative to HXB2
              try:
                with open(AlignmentFile, 'r') as file:
                  HXB2Alignment = AlignIO.read(file, 'fasta')
                  sequence_pos = 0
                  for j in range(len(HXB2Alignment[0].seq)):
                    if HXB2Alignment[0][j] != '-':
                      sequence_pos += 1
                    if sequence_pos == genome_indel_pos:
                      hxb2_alignment_pos = j
                      hxb2_indel_pos = len(ungap(HXB2Alignment[1].seq[0:int(hxb2_alignment_pos) + 1]))
                      indels["Indel{}".format(i)].update({"HXB2Position": hxb2_indel_pos})
                      break
              except IOError:
                print(('Error: Could not open file {}'.format(AlignmentFile)))
                raise
              except ValueError:
                print(('Error: File {} is not in the expected format'.format(AlignmentFile)))
                raise
              
      # Write indels to CSV
      with open(OutputCSV, 'a', newline='') as csvfile:
        fieldnames = ['Sequence', 'Gene', 'IndelID', 'Type', 'GenePosition', 'GenomePosition', 'GlobalGenomeRefPosition', 'HXB2Position']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        if csvfile.tell() == 0:
          writer.writeheader()
        for indel_id, data in list(indels.items()):
          writer.writerow({'Sequence': data['Sequence'], 'Gene': data['Gene'], 'IndelID': indel_id, 'Type': data['Type'], 'GenePosition': data['GenePosition'], 'GenomePosition': data['GenomePosition'], 'GlobalGenomeRefPosition': data['GlobalGenomeRefPosition'], 'HXB2Position': data['HXB2Position']})

  except IOError:
    print(('Error: Could not open file {}'.format(InputFile)))
    raise
  except ValueError:
    print(('Error: File {} is not in the expected format'.format(InputFile)))
    raise
  return 0

# Call CategoriseIndels
if function_name == 'CategoriseIndels':
  CategoriseIndels(args.InputFile, args.OutputCSV, args.Gene, args.SequenceName, args.ReferenceName, args.GeneCoordInfo, args.GenomeFile, args.AlignmentFile, args.ExtractedGeneFile)
