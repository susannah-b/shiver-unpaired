#!/usr/bin/env python
from __future__ import print_function

## Author: Chris Wymant, c.wymant@imperial.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''This script, taking an sequence alignment in fasta
format as input, considers all possible pairs of sequences therein and
calculates the sizes and positions (with respect to the alignment) of indels
if those two sequences were aligned on their own, i.e. not counting any position
at which both sequences have a gap.'''

import argparse
import os
import sys
import itertools
from Bio import AlignIO
from Bio import Seq  
from Bio import SeqIO  
from collections import Counter
import numpy as np


# Define a function to check files exist, as a type for the argparse.
def File(MyFile):
  if not os.path.isfile(MyFile):
    raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
  return MyFile

# Set up the arguments for this script
parser = argparse.ArgumentParser(description=ExplanatoryMessage)
parser.add_argument('alignment', type=File)
parser.add_argument('OutputFileBasename')
args = parser.parse_args()


# Read in the alignment
try:
  alignment = AlignIO.read(args.alignment, "fasta")
except:
  print('Problem reading', args.alignment + ':', file=sys.stderr)
  raise
AlignmentLength = alignment.get_alignment_length()

# Check there are 2+ sequences
NumSeqs = len(alignment)
if NumSeqs < 2:
  print('Need at least two sequences. Quitting.', file=sys.stderr)
  exit(1)
NumComparisons = (NumSeqs * (NumSeqs - 1)) / 2

# Our representation of each seq will be an array of bools: True for a gap char,
# False otherwise
seqs = np.array([np.array([base == '-' for base in seq.seq]) for seq in \
alignment])

def GetSeqStartAndEndPos(seq):
  '''Get the position of the first and last non-gap character in the seq.'''
  FirstBasePos = 0
  try:
    while seq[FirstBasePos]:
      FirstBasePos += 1
  except IndexError:
    print('Encountered pure-gap sequence. Quitting', file=sys.stderr)
    quit(1)
  LastBasePos = len(seq) - 1
  while seq[LastBasePos]:
    LastBasePos -= 1
  return FirstBasePos, LastBasePos

StartsAndEnds = []
for seq in seqs:
  start, end = GetSeqStartAndEndPos(seq)
  StartsAndEnds.append((start, end))

def ProcessRangeOfSeqs(SeqNumbers):
  '''Counts the del sizes and start positions for all pairwise comparisons
  between the specified seq numbers and greater seq numbers.'''

  DelSizeCounts = Counter()
  DelPositionCounts = Counter()
  count = 0

  for i in SeqNumbers:
    seq1 = seqs[i]
    start1, end1 = StartsAndEnds[i]

    for j, seq2 in enumerate(seqs[i+1:]):
      j += i+1
      start2, end2 = StartsAndEnds[j]

      count +=1
      if count % 2000 == 0:
        print('Now working on pair', count, 'of', NumComparisons)

      # Restrict our focus to the bit of the alignment after the start of both
      # seqs and before the end of either seq.
      start = max(start1, start2)
      end = min(end1, end2)
      if start > end:
        print("Encountered two seqs that don't overlap. Quitting",
        file=sys.stderr)
        exit(1)
      seq1sub = seq1[start : end + 1]
      seq2sub = seq2[start : end + 1]

      for pos, (gap1, gap2) in enumerate(itertools.izip(seq1sub, seq2sub)):

        # The first position in the pairwise alignment.
        if not pos:
          assert not (gap1 and gap2), "Malfunction: both sequences start with"+\
          " gaps after we've tried to start where the later one has its " + \
          "first base."
          if gap1:
            EffDelSize1 = 1
            DelStartPos1 = pos + start
          if gap2:
            EffDelSize2 = 1
            DelStartPos2 = pos + start
          LastGap1 = gap1
          LastGap2 = gap2
          continue

        if gap1:
          if LastGap1:
            # We're still inside a gap. If EffDelSize == 0, the other seq has
            # also had gaps for the whole deletion so far. If that's the case
            # and the other seq has a base now, the deletion start pos is here.
            # Also, increment EffDelSize by 1 if the other seq has a base
            if not EffDelSize1 and not gap2:
              DelStartPos1 = pos + start
            EffDelSize1 += not gap2
          else:
            # Starting a new gap!
            EffDelSize1 = int(not gap2)
            if EffDelSize1:
              DelStartPos1 = pos + start
        elif LastGap1:
          # this is the end of a gap
          if EffDelSize1:
            # It was a non-trivial gap - the other seq had at least one base
            # inside it.
            DelPositionCounts[DelStartPos1] += 1
            DelSizeCounts[EffDelSize1] += 1

        if gap2:
          if LastGap2:
            if not EffDelSize2 and not gap1:
              DelStartPos2 = pos + start
            EffDelSize2 += not gap1
          else:
            EffDelSize2 = int(not gap1)
            if EffDelSize2:
              DelStartPos2 = pos + start
        elif LastGap2:
          if EffDelSize2:
            DelPositionCounts[DelStartPos2] += 1
            DelSizeCounts[EffDelSize2] += 1

        LastGap1 = gap1
        LastGap2 = gap2

  return DelSizeCounts, DelPositionCounts
  
# Do all the pairwise comparisons!
DelSizeCounts, DelPositionCounts = ProcessRangeOfSeqs(range(NumSeqs))

def FillInCounterBlanksAndRescale(counter):
  '''Rescale all values by n(n-1)/2 and supply zero for missing values.'''
  MaxValSize = max(counter.keys())
  for val in range(1, MaxValSize + 1):
    if val in counter:
      counter[val] = float(counter[val]) / NumComparisons
    else:
      counter[val] = 0

FillInCounterBlanksAndRescale(DelPositionCounts)
FillInCounterBlanksAndRescale(DelSizeCounts)

# Write output  
with open(args.OutputFileBasename + '.csv', 'w') as f:
  f.write('Indel size (bp),Mean number of indels of that size per compared pair of sequences\n')
  for DelSize, DelSizeCount in DelSizeCounts.items():
    f.write(str(DelSize) + ',' + str(DelSizeCount) + '\n')
with open(args.OutputFileBasename + '_positions.csv', 'w') as f:
  f.write('Indel pos in alignment,Count\n')
  for DelPosition, DelPositionCount in DelPositionCounts.items():
    # Shift positions from 0-based to 1-based!
    f.write(str(DelPosition + 1) + ',' + str(DelPositionCount) + '\n')

