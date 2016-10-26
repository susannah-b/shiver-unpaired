#!/usr/bin/env python
from __future__ import print_function

## Author: Chris Wymant, c.wymant@imperial.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''Used with two base frequency files of the format
produced by shiver, each one resulting from mapping to a different reference,
and an alignment which contains those two references, this script aligns the two
coverages to each other based on the alignment. i.e., at each position in the
alignment, we report the coverage (number of mapped reads) with respect to each
reference.'''

import argparse
import os
import sys
from Bio import AlignIO
import itertools

# Define a function to check files exist, as a type for the argparse.
def File(MyFile):
  if not os.path.isfile(MyFile):
    raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
  return MyFile

# Set up the arguments for this script
parser = argparse.ArgumentParser(description=ExplanatoryMessage)
parser.add_argument('alignment', type=File, help='An alignment that contains '+\
'both references.')
parser.add_argument('ref1name')
parser.add_argument('ref2name')
parser.add_argument('BaseFreqsFile1', type=File)
parser.add_argument('BaseFreqsFile2', type=File)
args = parser.parse_args()

# Check for unique ref names
if args.ref1name == args.ref2name:
  print('You cannot specify the same reference twice. Quitting.', \
  file=sys.stderr)
  exit(1)

# Find the consensus and its ref
alignment = AlignIO.read(args.alignment, "fasta")
ref1seq = None
ref2seq = None
for seq in alignment:
  if seq.id == args.ref1name:
    if ref1seq != None:
      print('Found', args.ref1name, 'twice in', args.alignment + '. Quitting.',\
      file=sys.stderr)
      quit(1)
    ref1seq = str(seq.seq)
  if seq.id == args.ref2name:
    if ref2seq != None:
      print('Found', args.ref2name, 'twice in', args.alignment + '. Quitting.',\
      file=sys.stderr)
      quit(1)
    ref2seq = str(seq.seq)

# Check the refs were found
if ref1seq == None:
  print(args.ref1name, 'not found in', args.alignment + \
  '. Quitting.', file=sys.stderr)
  exit(1)
if ref2seq == None:
  print(args.ref2name, 'not found in', args.alignment + \
  '. Quitting.', file=sys.stderr)
  exit(1)

def GetCoverages(BaseFreqsFile, RefSeq):
  '''Converts coverages with respect to the ref to coverages with respect to the
  alignment.'''

  GaplessRefSeq = RefSeq.replace('-','').upper()
  CoveragesInRef = []
  with open(BaseFreqsFile, 'r') as f:
    LastRefPos = 0
    for LineNumMin1, line in enumerate(f):
      if LineNumMin1 == 0:
        continue
      fields = line.split(',')
      RefPos = fields[0]
      if RefPos == 'NA':  
        continue
      RefPos = int(RefPos)
      if RefPos != LastRefPos + 1:
        print('Error: the reference position jumped to', RefPos, 'on line', \
        LineNumMin1+1, 'in', BaseFreqsFile + '. Quitting', file=sys.stderr)
        quit(1)
      RefBase = fields[1]
      if RefBase.upper() != GaplessRefSeq[RefPos-1]:
        print('Error: mismatch in reference base at position', str(RefPos) + \
        ':', RefBase, 'in', BaseFreqsFile, 'but', GaplessRefSeq[RefPos-1], \
        'in', args.alignment + '. Quitting.', file=sys.stderr)
        quit(1)
      coverage = sum(map(int, fields[2:]))
      CoveragesInRef.append(coverage)
      LastRefPos += 1
  RefLength = len(GaplessRefSeq)
  assert len(CoveragesInRef) == RefLength
  CoveragesInAlignment = []
  AlnPosToRefPos = []
  LastPosCoverage = 0
  ZeroBasedPosInRef = -1
  for base in RefSeq:
    if ZeroBasedPosInRef == RefLength - 1:
      # We're after the end of the ref in the alignment
      coverage = 0
    elif base == '-':
      # Zero if the ref hasn't started yet, else equal to its last value before
      # the insertion:
      coverage = LastPosCoverage
    else:
      ZeroBasedPosInRef += 1
      coverage = CoveragesInRef[ZeroBasedPosInRef]
    CoveragesInAlignment.append(coverage)
    LastPosCoverage = coverage
    if base == '-':
      AlnPosToRefPos.append('-')
    else:
      AlnPosToRefPos.append(ZeroBasedPosInRef+1)
  return CoveragesInAlignment, AlnPosToRefPos

ref1coverages, ref1PosConversions = GetCoverages(args.BaseFreqsFile1, ref1seq)
ref2coverages, ref2PosConversions = GetCoverages(args.BaseFreqsFile2, ref2seq)

print('Alignment position,Position in ' + args.ref1name + ',Position in ' + args.ref2name + ',Coverage for ' + args.ref1name + ',Coverage for ' + \
args.ref2name)
for PosMin1, (ref1cov, ref2cov) in enumerate(itertools.izip(ref1coverages, \
ref2coverages)):
  PosInRef1 = ref1PosConversions[PosMin1]
  PosInRef2 = ref2PosConversions[PosMin1]
  print(PosMin1+1, PosInRef1, PosInRef2, ref1cov, ref2cov, sep=',')
