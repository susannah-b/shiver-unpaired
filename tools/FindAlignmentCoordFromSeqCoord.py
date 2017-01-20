#!/usr/bin/env python
from __future__ import print_function

## Author: Chris Wymant, c.wymant@imperial.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''This alignment translates positions with respect to a
named sequence to positions with respect to an alignment containing that
sequence. e.g. if your named seq started like this in the alignment: --a-g...
and you specified position 2, this script will tell you position 5.'''

import argparse
import os
import sys
from Bio import AlignIO

# Define a function to check files exist, as a type for the argparse.
def File(MyFile):
  if not os.path.isfile(MyFile):
    raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
  return MyFile

# Set up the arguments for this script
ExplanatoryMessage = ExplanatoryMessage.replace('\n', ' ').replace('  ', ' ')
parser = argparse.ArgumentParser(description=ExplanatoryMessage)
parser.add_argument('alignment', type=File)
parser.add_argument('SeqName')
parser.add_argument('PositionsInSeq', type=int, nargs='+',
help='One or more positive integers.')
args = parser.parse_args()

# Check all positions positive
if min(args.PositionsInSeq) < 1:
  print('All positions should be greater than zero. Quitting.', file=sys.stderr)
  exit(1)

alignment = AlignIO.read(args.alignment, "fasta")
ChosenSeq = None
for seq in alignment:
  if seq.id == args.SeqName:
    ChosenSeq = str(seq.seq)
    break
if ChosenSeq == None:
  print('Did not find', args.SeqName, 'in', args.alignment + '. Quitting.', 
  file=sys.stderr)
  exit(1)

# Count missing coverage as a gap
ChosenSeq = ChosenSeq.replace('?', '-')

# Check no positions are after the end of the seq
if max(args.PositionsInSeq) > len(ChosenSeq) - ChosenSeq.count('-'):
  print('You specified position', max(args.PositionsInSeq),
  'but', args.SeqName, 'is only', len(ChosenSeq) - ChosenSeq.count('-'),
  'bases long. Quitting.', file=sys.stderr)
  exit(1)

# Stolen from phyloscanner
def TranslateSeqCoordsToAlnCoords(seq, coords):
  '''Takes a sequence that contains gaps (in general), and a set of coordinates
  specified with a respect to that sequence without gaps. The coordinates are
  translated to their positions in the gappy version of the sequence.
  e.g. called with the arguments "-a--cg-t-" and [1,2,3], we return [2,5,6].
  '''
  TranslatedCoords = [-1 for coord in coords]
  PositionInSeq = 0
  for GappyPostitionMin1,base in enumerate(seq):
    if base != '-':
      PositionInSeq += 1
      for i,coord in enumerate(coords):
        if coord == PositionInSeq:
          TranslatedCoords[i] = GappyPostitionMin1+1
      if not -1 in TranslatedCoords:
        break
  assert not -1 in TranslatedCoords
  assert len(TranslatedCoords) == len(coords)
  return TranslatedCoords

TranslatedCoords = TranslateSeqCoordsToAlnCoords(ChosenSeq, args.PositionsInSeq)
print(' '.join(map(str, TranslatedCoords)))
