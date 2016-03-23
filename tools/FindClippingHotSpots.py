#!/usr/bin/env python
from __future__ import print_function

## Author: Chris Wymant, c.wymant@imperial.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''TODO:
'''

import os
import collections
import sys
import argparse
import pysam

# Define a function to check files exist, as a type for the argparse.
def File(MyFile):
  if not os.path.isfile(MyFile):
    raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
  return MyFile

# Set up the arguments for this script
ExplanatoryMessage = ExplanatoryMessage.replace('\n', ' ').replace('  ', ' ')
parser = argparse.ArgumentParser(description=ExplanatoryMessage)
parser.add_argument('BamFile', type=File)
args = parser.parse_args()

BamFile = pysam.AlignmentFile(args.BamFile, "rb")

# Find the reference in the bam file; there should only be one.
AllReferences = BamFile.references
if len(AllReferences) != 1:
  print('Expected exactly one reference in', BamFileName+'; found',\
  str(len(AllReferences))+'.\nQuitting.', file=sys.stderr)
  exit(1)
RefName = AllReferences[0]

# Get the length of the reference.
AllReferenceLengths = BamFile.lengths
if len(AllReferenceLengths) != 1:
  print('Pysam error: found one reference but', len(AllReferenceLengths), \
  'reference lengths.\nQuitting.', file=sys.stderr)
  exit(1)
RefLength = AllReferenceLengths[0]

ClipPositions = []
NumbersOfSpanningReads = [0 for i in range(RefLength)]
NumReads = 0

for read in BamFile.fetch(RefName):

  positions = read.get_reference_positions(full_length=True)

  # Shouldn't happen, but skip reads mapped to no position
  if not any(pos != None for pos in positions):
    continue

  NumReads += 1

  # If the left edge is clipped, find where.
  LeftMostMappedBase = 0
  while positions[LeftMostMappedBase] == None:
    LeftMostMappedBase += 1
  if LeftMostMappedBase > 0:
    ClipPositions.append(positions[LeftMostMappedBase])

  # If the right edge is clipped, find where.
  RightMostMappedBase = len(positions)-1
  while positions[RightMostMappedBase] == None:
    RightMostMappedBase -= 1
  if RightMostMappedBase < len(positions)-1:
    ClipPositions.append(positions[RightMostMappedBase]+1)

  # Update the counts of reads spanning each position
  for pos in range(positions[LeftMostMappedBase]+1, \
  positions[RightMostMappedBase]+1):
    NumbersOfSpanningReads[pos] += 1

if NumReads == 0:
  print('No mapped reads found in', args.BamFile, '. Quitting.', \
  file=sys.stderr)
  exit(1)

# Count how many times each position was recorded
ClipPositionCounts = collections.Counter(ClipPositions)

# Print the output
print('Reference position, Number of reads clipped, Percentage of spanning',\
'reads clipped')
for pos, count in sorted(ClipPositionCounts.items(), key=lambda x:x[1], \
reverse=True):
  NumberOfSpanningReads = NumbersOfSpanningReads[pos]
  PercentageClipped = 100 * float(count) / (count + NumberOfSpanningReads)
  print(pos+1, count, '%.3f' % PercentageClipped, sep=',')

