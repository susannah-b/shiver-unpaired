#!/usr/bin/python
from __future__ import print_function

## Author: Chris Wymant, chris.wymant@bdi.ox.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''TODO'''

import argparse
import os
import sys
import itertools
from Bio import AlignIO
from Bio import Seq  
from Bio import SeqIO  
import collections
from re import sub
import numpy as np
import matplotlib.pyplot as plt

# Define a function to check files exist, as a type for the argparse.
def File(_file):
  if not os.path.isfile(_file):
    raise argparse.ArgumentTypeError(_file +' does not exist or is not a file.')
  return _file


# Set up the arguments for this script
parser = argparse.ArgumentParser(description=ExplanatoryMessage)
parser.add_argument('alignment', type=File)
parser.add_argument('patient_based_blacklist', type=File)
parser.add_argument('seq_based_blacklist', type=File)
parser.add_argument('amplicon_regions_file', type=File)
parser.add_argument('output_file')
parser.add_argument('-V', '--verbose', action='store_true',
help="Print a line for each blacklisting action.")
parser.add_argument('-SA', '--split-amplicons', action='store_true')
parser.add_argument('--print_new_seq_blacklist', help='''Use this to specify the
name of file in which we'll store additions to the seq_based_blacklist (due to
our analysis here of completeness in each region).''')
args = parser.parse_args()

# Read the alignment
try:
  alignment = AlignIO.read(args.alignment, "fasta")
except:
  print('Problem reading', args.alignment + ':', file=sys.stderr)
  raise
alignment_length = alignment.get_alignment_length()

def get_beehive_id(seq_id):
  "Return whatever's before the first underscore if there is one."
  return seq_id.split('_', 1)[0]

# Read in the seqs, replace "?" and lower case letters by "N", check for
# duplicate seqs and unexpected bases.
seq_dict = collections.OrderedDict()
for seq in alignment:
  id_ = seq.id
  id_ = id_.split('_consensus', 1)[0]
  id_ = id_.split('_MinCov', 1)[0]
  if id_ in seq_dict:
    print('Encountered seq', id_, 'a second time in',
    args.alignment + '. Quitting.', file=sys.stderr)
    exit(1)
  SeqAsStr = str(seq.seq)
  SeqAsStr = sub("[a-z]|\?", "N", SeqAsStr)
  assert len(SeqAsStr) == alignment_length
  if any(not base in "ACGTN-" for base in SeqAsStr):
    print('Seq', id_, 'in', args.alignment, 'contains a base other than A, C,',
    'G, T, N or -; unexpected. Have you run',
    'shiver/tools/EstimateAmbiguousBases.py on your alignment first? Quitting.',
    file=sys.stderr)
    exit(1)
  if SeqAsStr.count("N") + SeqAsStr.count("-")  == alignment_length:
    if args.verbose:
      print("Skipping sequence", id_, "which is wholly undetermined after",
      "removing lower-case bases.")
  else:
    seq_dict[id_] = SeqAsStr

all_seq_ids = set(seq_dict.keys())
all_patients_with_seqs = set([get_beehive_id(seq_id) for seq_id in \
seq_dict])

# Read and check the amplicon regions. 
regions_dict = collections.OrderedDict()
regions_dict_not_empty = False
with open(args.amplicon_regions_file, 'r') as f:
  for line in f:
    fields = line.split(',')
    if len(fields) != 3:
      print('Encountered', len(fields), 'fields in', \
      args.amplicon_regions_file + '. Expected 3. Quitting.', file=sys.stderr)
      exit(1)
    region, start, end = fields
    start, end = int(start), int(end)
    if start > end:
      print('Start greater than end for region', region, 'in', 
      args.amplicon_regions_file + '. Quitting.', file=sys.stderr)
      exit(1)
    if region in regions_dict:
      print('Encountered', region, 'a second time in',
      args.amplicon_regions_file + '. Quitting.', file=sys.stderr)
      exit(1)
    if regions_dict_not_empty:
      previous_region = next(reversed(regions_dict))
      previous_start, previous_end = regions_dict[previous_region]
      if start != previous_end + 1:
        print('Region', region, 'starts at', start, 'which is not 1 more than',
        'the end of the previous region (' + str(previous_end), 'for',
        previous_region + '). Quitting.', file=sys.stderr)
        exit(1)
    elif start != 1:
      print('The first region in', args.amplicon_regions_file,
      'does not start at 1. Quitting.', file=sys.stderr)
      exit(1)
    regions_dict[region] = (start, end)
    regions_dict_not_empty = True
last_region = next(reversed(regions_dict))
last_start, last_end = regions_dict[last_region]
if last_end != alignment_length:
  print('The last region in', args.amplicon_regions_file,
  'does not end at the alignment length (' + str(alignment_length) + \
  '). Quitting.', file=sys.stderr)
  exit(1)
regions = regions_dict.keys()
num_regions = len(regions)

# Read the seq-based blacklist.
seq_blacklist_dict = {}
expected_header_line = 'BAM,keep.at.all,' + ','.join(regions) + ',origin'
num_fields = expected_header_line.count(',') + 1
with open(args.seq_based_blacklist, 'r') as f:
  for line_num_min_1, line in enumerate(f):

    # Check for the expected header
    if line_num_min_1 == 0:
      line = line.strip()
      if line != expected_header_line:
        print('Unexpected header line\n' + line + '\nfor',
        args.seq_based_blacklist + '; expected\n' + \
        expected_header_line + '\nQuitting.', file=sys.stderr)
        exit(1)
      continue

    # Check for the right number of fields
    fields = line.split(',')
    if len(fields) != num_fields:
      print('Unexpected number of fields', len(fields), 'on line',
      line_num_min_1 + 1, 'in', args.seq_based_blacklist + '; expected',
      str(num_fields) + '. Quitting.', file=sys.stderr)
      exit(1)

    # Warn about, and skip, blacklisted sequences that do not appear in the
    # global alignment.
    id_ = fields[0]
    if not id_ in all_seq_ids:
      print("Skipping sequence", id_, "from",
      args.seq_based_blacklist, "which does not appear in", args.alignment + \
      ".", file=sys.stderr)
      continue

    # Coerce string bools to bools, and record them.
    values = fields[1:-1]
    if any(value != "TRUE" and value != "FALSE" for value in values):
      print('Unexpected value not equal to TRUE or FALSE on line',
      line_num_min_1 + 1, 'in', args.seq_based_blacklist + '. Quitting.',
      file=sys.stderr)
      exit(1)
    values = np.array([True if value == "TRUE" else False for value in values],
    dtype=bool)

    # If we've seen this seq in the blacklist already, we should blacklist each
    # region if either entry says so, i.e. only keep it if they agree that we
    # should.
    if id_ in seq_blacklist_dict:
      previous_values = seq_blacklist_dict[id_]
      for i, new_value in enumerate(values):
        previous_value = previous_values[i]
        if new_value and previous_value:
          values[i] = True
        else:
          values[i] = False

    seq_blacklist_dict[id_] = values

blacklisted_patients = set([])
with open(args.patient_based_blacklist, 'r') as f:
  for line in f:
    patient = line.strip()
    if not patient in all_patients_with_seqs:
      print("Skipping blacklisted patient", patient, "from",
      args.patient_based_blacklist, "who does not appear in", args.alignment + \
      ".", file=sys.stderr)
    else:
      blacklisted_patients.add(patient)

#completeness_percents = collections.defaultdict(list)
max_missingness = 0.2
extra_blacklisted_seq_ids = set([])
for seq_id in all_seq_ids:

  # Add our own blacklisting of regions, based on the fraction that's not "N",
  # for regions that have not already been blacklisted.
  for region_num in xrange(num_regions):
    seq_in_blacklist = seq_id in seq_blacklist_dict
    if seq_in_blacklist and not seq_blacklist_dict[seq_id][region_num + 1]:
      continue
    region = regions[region_num]
    start, end = regions_dict[region]
    region_length = end - start + 1
    seq_here = seq_dict[seq_id][start - 1: end]
    missingness = float(seq_here.count("N")) / region_length
    if missingness >= max_missingness:
      if not seq_in_blacklist:
        seq_blacklist_dict[seq_id] = [True] * (num_regions + 1)
      seq_blacklist_dict[seq_id][region_num + 1] = False
      extra_blacklisted_seq_ids.add(seq_id)
    #completeness_percents[region].append(completeness_percent)

  # Delete every seq from a blacklisted patient
  beehive_id = get_beehive_id(seq_id)
  if beehive_id in blacklisted_patients:
    if args.verbose:
      print("Sequence ", seq_id, ": removing whole sequence as patient ",
      beehive_id, " was blacklisted.", sep='')
    del seq_dict[seq_id]
    continue

  if seq_id in seq_blacklist_dict:

    # The seq_blacklist_values are bools saying keep this seq at all, then keep
    # each region.
    seq_blacklist_values = seq_blacklist_dict[seq_id]
    keep_at_all = seq_blacklist_values[0]
    if not keep_at_all:
      if args.verbose:
        print("Sequence ", seq_id,
        ": removing whole sequence as it was blacklisted.", sep='')
      del seq_dict[seq_id]
      continue

    for region_num in xrange(num_regions):

      # Mask blacklisted regions by "N".
      # Remember start and end are 1-based coords.
      keep_region = seq_blacklist_values[region_num + 1]
      if not keep_region:
        region = regions[region_num]
        start, end = regions_dict[region]
        if args.verbose:
          print("Sequence ", seq_id, ": masking region ", region + \
          " (positions ", start, "-", end, ").", sep='')
        seq_dict[seq_id] = seq_dict[seq_id][:start - 1] + \
        "N" * (end - start + 1) + seq_dict[seq_id][end:]

if args.print_new_seq_blacklist != None:
  with open(args.print_new_seq_blacklist, 'w') as f:
    f.write(expected_header_line.rsplit(',', 1)[0] + "\n")
    for seq_id in sorted(extra_blacklisted_seq_ids):
      seq_blacklist_values = seq_blacklist_dict[seq_id]
      f.write(seq_id + "," + ",".join(map(str, seq_blacklist_values)).upper() \
      + "\n")

#for region, completeness_percents_here in completeness_percents.items():
#  fig, ax = plt.subplots()
#  ax.set_yscale("log", nonposy='clip')
#  ax.set_ylim(bottom=0.5, top=len(alignment))
#  plt.xlabel("Percent of alignment positions that are not 'N'", fontsize=15)
#  plt.ylabel("Number of sequences", fontsize=15)
#  plt.hist(completeness_percents_here, bins=100)
#  plt.savefig("AmpliconCompletenessHisto_" + region + ".pdf")
#  plt.clf()
#exit(0)

# Form a dict for patients with multiple sequences, listing the sequence IDs
per_patient_seq_counts = \
collections.Counter(get_beehive_id(seq_id) for seq_id in seq_dict)
multiply_seqd_patients = set(patient for patient, count in \
per_patient_seq_counts.items() if count > 1)
singly_seqd_patients = set(patient for patient, count in \
per_patient_seq_counts.items() if count == 1)
seq_ids_for_multiply_seqd_patients = collections.defaultdict(list)
for seq_id in seq_dict:
  beehive_id = get_beehive_id(seq_id)
  if beehive_id in multiply_seqd_patients:
    seq_ids_for_multiply_seqd_patients[beehive_id].append(seq_id)


# Merge multiple seqs per patient.
for multiply_seqd_patient, seq_ids in \
seq_ids_for_multiply_seqd_patients.items():
  seqs = [seq_dict[seq_id] for seq_id in seq_ids]
  num_seqs = len(seqs)
  seqs_by_length = sorted(seqs, key=lambda seq : alignment_length - \
  seq.count("N") - seq.count("-"), reverse=True)
  best_seq = ""
  for pos in xrange(alignment_length):

    # Try each seq in order from best to worst until one of them has something
    # other than an N:
    best_base = "N"
    for i in xrange(num_seqs):
      base = seqs_by_length[i][pos]
      if base != "N":
        best_base = base
        break
    best_seq += best_base

  for seq_id in seq_ids:
    del seq_dict[seq_id]
  seq_dict[multiply_seqd_patient] = best_seq

# Multiply sequenced patients now have only one seq, labelled by the patient;
# ensure that for other patients their only seq is also labelled by the patient.
seq_ids_to_rename = \
[seq_id for seq_id in seq_dict if get_beehive_id(seq_id) != seq_id]
for seq_id in seq_ids_to_rename:
  seq_dict[get_beehive_id(seq_id)] = seq_dict.pop(seq_id)


def PropagateNoCoverageChar(seq, LeftToRightDone=False):
  '''Iteratively replace gaps that border N by N.

  Where N neighbours a gap, propagate N outwards until it touches bases on both
  sides (because deletions should only be called when the bases on either side
  are known). e.g.
  ACTG---N---ACTG
  becomes
  ACTGNNNNNNNACTG'''
  
  if LeftToRightDone:
    seq = seq[::-1]
  BaseToLeftIsNoCoverage = False
  ResultingSeq = ''
  for base in seq:
    if base == 'N':
      BaseToLeftIsNoCoverage = True
      ResultingSeq += 'N'
    elif base == '-':
      if BaseToLeftIsNoCoverage:
        ResultingSeq += 'N'
      else:
        ResultingSeq += '-'
    else:
      BaseToLeftIsNoCoverage = False
      ResultingSeq += base
  if LeftToRightDone:
    ResultingSeq = ResultingSeq[::-1]
  else:
    ResultingSeq = PropagateNoCoverageChar(ResultingSeq, True)
  assert not "N-" in ResultingSeq and not "-N" in ResultingSeq
  return ResultingSeq

# Iteratively replace gaps that border N by N.
for seq_id in seq_dict:
  seq_dict[seq_id] = PropagateNoCoverageChar(seq_dict[seq_id])

# Sort seqs by name
sorted_seqs = sorted(seq_dict.items(), key=lambda x:x[0])

if args.split_amplicons:

  # Set up file names for the per-region output files. Use an extension if the
  # user has specified one (splitting into two pieces at the right-most dot).
  pieces = args.output_file.rsplit(".", 1)
  if len(pieces) == 2:
    per_region_output_file_dict = {region:pieces[0] + "_" + region + "." + \
    pieces[1] for region in regions}
  else:
    per_region_output_file_dict = {region:args.output_file + "_" + region + \
    ".fasta" for region in regions}

  for region, (start, end) in regions_dict.items():
    OutSeqs = []
    region_length = end - start + 1
    max_unknown_chars = region_length / 2 # int division, OK
    for seq_id, seq in sorted_seqs:
      seq_here = seq[start - 1: end]
      if seq_here.count("N") <= max_unknown_chars:
        OutSeqs.append(SeqIO.SeqRecord(Seq.Seq(seq_here), id=seq_id,
        description=''))
    SeqIO.write(OutSeqs, per_region_output_file_dict[region], 'fasta')    

else:
  OutSeqs = []
  for seq_id, seq in sorted_seqs:
    if seq.count("N") == alignment_length:
      print("Skipping sequence", seq_id, "which is wholly undetermined after",
      "blacklisting.")
    else:
      OutSeqs.append(SeqIO.SeqRecord(Seq.Seq(seq), id=seq_id, description=''))
  SeqIO.write(OutSeqs, args.output_file, 'fasta')
        

