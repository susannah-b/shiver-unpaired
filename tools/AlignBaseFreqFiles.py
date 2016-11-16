#!/usr/bin/env python
from __future__ import print_function

## Author: Chris Wymant, c.wymant@imperial.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''Used with two base frequency files of the format
produced by shiver, each one resulting from mapping to a different reference,
and an alignment which contains those two references, this script aligns the two
frequencies to each other based on the alignment. i.e., at each position in the
alignment, we report the frequencies with respect to each reference. Output is
printed to stdout suitable for redirection to a csv file.'''

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
parser.add_argument('-C', '--coverage-only', action='store_true', help='''At
each position in the alignment, instead of printing the base frequencies for
each reference, print their coverage (the total number of reads mapped there,
i.e. the sum of the base counts).''')
parser.add_argument('-CF', '--compare', action='store_true', help='''At each
position in the alignment, instead of printing the base frequencies for both
references, print a score in the range [0,1] to summarise how similar the two
sets of frequencies are to each other. That score is calculated as follows. The
count for "N" is ignored. For each base "A", "C", "G", "T", "-", we calculate
the absolute difference between the fraction of all reads each reference has
with that base; we then sum these values, divide by 2, and subtract the
result from 1. This means a score of 1 is obtained if and only if the two sets
of frequencies agree perfectly, a score of zero is obtained if and only if there
are no bases in common, and all other situations give an intermediate value. A
score of NA is given at positions where one or both references have no mapped
reads, or if one reference has an insertion with respect to the other.''')
parser.add_argument('-CS', '--compare-simple', action='store_true', help='''Like
--compare, except that the following simple similarity metric is used: 0 for
disagreement on which of the five bases is most common, 1 for agreement.''')
args = parser.parse_args()

# Check for unique ref names
if args.ref1name == args.ref2name:
  print('You cannot specify the same reference twice. Quitting.',
  file=sys.stderr)
  exit(1)

# No commas in reference names
if ',' in args.ref1name or ',' in args.ref2name:
  print('Reference names may not contain commas. Rename the sequences in',
  args.alignment, 'and try again. Quitting.', file=sys.stderr)
  exit(1)

# Can't use both --compare and --compare-simple
if args.compare and args.compare_simple:
  print('You cannot use both --compare and --compare-simple. Quitting.',
  file=sys.stderr)
  exit(1)

# Find the consensus and its ref
alignment = AlignIO.read(args.alignment, "fasta")
ref1seq = None
ref2seq = None
for seq in alignment:
  if seq.id == args.ref1name:
    if ref1seq != None:
      print('Found', args.ref1name, 'twice in', args.alignment + '. Quitting.',
      file=sys.stderr)
      quit(1)
    ref1seq = str(seq.seq)
  if seq.id == args.ref2name:
    if ref2seq != None:
      print('Found', args.ref2name, 'twice in', args.alignment + '. Quitting.',
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

def GetFreqs(BaseFreqsFile, RefSeq):
  '''Converts freqs with respect to the ref to freqs with respect to the
  alignment.'''

  GaplessRefSeq = RefSeq.replace('-','').upper()
  FreqsInRef = []
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
        print('Error: the reference position jumped to', RefPos, 'on line',
        LineNumMin1+1, 'in', BaseFreqsFile + '. Quitting', file=sys.stderr)
        quit(1)
      RefBase = fields[1]
      if RefBase.upper() != GaplessRefSeq[RefPos-1]:
        print('Error: mismatch in reference base at position', str(RefPos) + \
        ':', RefBase, 'in', BaseFreqsFile, 'but', GaplessRefSeq[RefPos-1],
        'in', args.alignment + '. Quitting.', file=sys.stderr)
        quit(1)
      freqs = map(int, fields[2:])
      assert len(freqs) == 6
      FreqsInRef.append(freqs)
      LastRefPos += 1
  RefLength = len(GaplessRefSeq)
  assert len(FreqsInRef) == RefLength
  FreqsInAlignment = []
  AlnPosToRefPos = []
  LastPosFreqs = [0,0,0,0,0,0]
  ZeroBasedPosInRef = -1
  for base in RefSeq:
    if ZeroBasedPosInRef == RefLength - 1:
      # We're after the end of the ref in the alignment
      freqs = [0,0,0,0,0,0]
    elif base == '-':
      # If the ref hasn't started yet, all base counts equal zero; if this is a
      # gap inside the ref, reproduce the base counts from the last position.
      # This is just to facilitate the coverage calculation, for which we take 
      # the coverage (sum of counts) here to be equal to its last value before
      # the deletion, but we won't report the breakdown into different bases.
      freqs = LastPosFreqs
    else:
      ZeroBasedPosInRef += 1
      freqs = FreqsInRef[ZeroBasedPosInRef]
    FreqsInAlignment.append(freqs)
    LastPosFreqs = freqs
    if base == '-':
      AlnPosToRefPos.append('-')
    else:
      AlnPosToRefPos.append(ZeroBasedPosInRef+1)
  return FreqsInAlignment, AlnPosToRefPos

ref1freqs, ref1PosConversions = GetFreqs(args.BaseFreqsFile1, ref1seq)
ref2freqs, ref2PosConversions = GetFreqs(args.BaseFreqsFile2, ref2seq)

# Set up the csv file header
outstring = 'Alignment position,Position in ' + args.ref1name + \
',Position in ' + args.ref2name
if args.coverage_only:
  outstring += ',coverage for ' + args.ref1name + ',coverage for ' + \
  args.ref2name
if args.compare:
  outstring += ',base frequency similarity score'
elif args.compare_simple:
  outstring += ',agreement on the most common base?'
elif not args.coverage_only:
  outstring += ',A count for ' + args.ref1name + \
  ',C count for ' + args.ref1name + ',G count for ' + args.ref1name + \
  ',T count for ' + args.ref1name + ',gap count for ' + args.ref1name + \
  ',N count for ' + args.ref1name + ',A count for ' + args.ref2name + \
  ',C count for ' + args.ref2name + ',G count for ' + args.ref2name + \
  ',T count for ' + args.ref2name + ',gap count for ' + args.ref2name + \
  ',N count for ' + args.ref2name

# Record each row of the csv file
for PosMin1, (ref1freqs, ref2freqs) in enumerate(itertools.izip(ref1freqs,
ref2freqs)):
  PosInRef1 = ref1PosConversions[PosMin1]
  PosInRef2 = ref2PosConversions[PosMin1]
  outstring += '\n' + str(PosMin1+1) + ',' + str(PosInRef1) + ',' + \
  str(PosInRef2) 

  if args.coverage_only:
    ref1cov = sum(ref1freqs)
    ref2cov = sum(ref2freqs)
    outstring += ',' + str(ref1cov) + ',' + str(ref2cov)

  if args.compare or args.compare_simple:
    ref1freqs = ref1freqs[:5]
    ref2freqs = ref2freqs[:5]  
    ref1cov = sum(ref1freqs)
    ref2cov = sum(ref2freqs)
    if ref1cov == 0 or ref2cov == 0 or PosInRef1 == '-' or PosInRef2 == '-':
      SimScore = 'NA'
    else:
      if args.compare_simple:
        MaxFreq1 = max(ref1freqs)
        BasesWithMaxCount1 = [i for i,count in enumerate(ref1freqs) \
        if count == MaxFreq1]
        MaxFreq2 = max(ref2freqs)
        BasesWithMaxCount2 = [i for i,count in enumerate(ref2freqs) \
        if count == MaxFreq2]
        if BasesWithMaxCount1 == BasesWithMaxCount2:
          SimScore = 1
        else:
          SimScore = 0
      else:
        SimScore = 0
        for i in range(5):
          SimScore += \
          abs(float(ref1freqs[i])/ref1cov - float(ref2freqs[i])/ref2cov)
        SimScore = 1 - SimScore/2
      SimScore = str(SimScore)
    outstring += ',' + SimScore

  elif not args.coverage_only:
    if PosInRef1 == '-':
      ref1freqs = ['-'] * 6
    if PosInRef2 == '-':
      ref2freqs = ['-'] * 6
    outstring += ',' + ','.join(map(str,ref1freqs)) + ',' + \
    ','.join(map(str,ref2freqs))

print(outstring)
