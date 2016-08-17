#!/usr/bin/env python2
from __future__ import print_function
#
## Author: Chris Wymant, c.wymant@imperial.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''
This script interprets a base frequency file of the format produced by
AnalysePileup.py and calls the consensus sequence, which is printed to stdout
suitable for redirection to a fasta-format file. The consensus is printed with
the reference used for mapping, as a pairwise alignment.
'''
##
################################################################################
# USER INPUT
# Gap characters in the base frequency file.
GapChar = '-'
################################################################################


# Import some modules we'll need.
import os.path
import sys
import itertools
from AuxiliaryFunctions import ReverseIUPACdict2, PropagateNoCoverageChar
import argparse
from Bio import SeqIO
from Bio import Seq


# Define a function to check files exist, as a type for the argparse.
def File(MyFile):
  if not os.path.isfile(MyFile):
    raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
  return MyFile

# Set up the arguments for this script
parser = argparse.ArgumentParser(description=ExplanatoryMessage)
parser.add_argument('BaseFreqFile', type=File)
parser.add_argument('MinCoverage', help='The minimum coverage (number of ' + \
'reads at a given position in the genome) before a base is called. Below ' + \
'this we call "?" instead of a base.', type=int)
parser.add_argument('MinCovForUpper', help='The minimum coverage before upper'+\
' case is used instead of lower case, to signal increased confidence.', \
type=int)
parser.add_argument('MinFracToCall', help='The minimum fraction of reads at a'+\
' position before we call that base (or those bases, when one base alone does'+\
' not reach that threshold fraction; e.g. say you have 60%% A, 30%% C and ' +\
'10%% G: if you set this fraction to 0.6 or lower we call an A, if you set ' +\
'it to 0.6-0.9 we call an M for "A or C", if you set it to 0.9-1 we call a ' + \
'V for "A, C or G".). Alternatively, if you choose a negative value, we '+\
'always call the single most common base regardless of its fraction, unless ' +\
'two or more bases are equally (most) common, then we call the ' + \
'ambiguity code for those bases.', type=float)
parser.add_argument('-C', '--consensus-seq-name', help='The name used for the'+\
' consensus in the fasta-format output (default: "consensus").', \
default='consensus')
parser.add_argument('-R', '--ref-seq-name', help='The name used for the '+\
'reference in the fasta-format output (default: "MappingReference").', \
default='MappingReference')
args = parser.parse_args()

BaseFreqFile = args.BaseFreqFile
MinCoverage = args.MinCoverage
MinCovForUpper = args.MinCovForUpper

# Check that MinCoverage and MinCovForUpper are positive integers, the 
# latter not smaller than the former.
if MinCoverage < 1:
  print('The specified MinumumCoverageToCallBase of', MinCoverage, \
  'is less than 1. Quitting.', file=sys.stderr)
  exit(1)
if MinCovForUpper < MinCoverage:
  print('The specified MinumumCoverageToUseUpperCase of', MinCoverage, \
  'is less than the specified MinumumCoverageToCallBase. Quitting.', \
  file=sys.stderr)
  exit(1)

# MinFracToCall should be <= 1 and != 0
if args.MinFracToCall > 1:
  print('MinFracToCall cannot be greater than 1. Quitting.', file=sys.stderr)
  exit(1)
FloatTolerance = 1e-5
if abs(args.MinFracToCall) < FloatTolerance:
  print('MinFracToCall should not equal zero. Quitting.', file=sys.stderr)
  exit(1)

CallMostCommon = args.MinFracToCall < 0

def CallAmbigBaseIfNeeded(bases, coverage):
  '''If several bases are supplied, calls an ambiguity code. Uses upper/lower
  case appropriately.'''

  bases = ''.join(sorted(bases))
  NumBases = len(bases)
  assert NumBases > 0, 'CallAmbigBaseIfNeeded function called with no bases.'
  if len(bases) == 1:
    BaseHere = bases
  else:
    # If a gap is one of the things most common at this position, call an 'N';
    # otherwise, find the ambiguity code for this set of bases.
    if GapChar in bases:
      BaseHere = 'N'
    else:  
      try:
        BaseHere = ReverseIUPACdict2[bases]
      except KeyError:
        print('Unexpected set of bases', bases, 'found in', BaseFreqFile, \
        ', not found amonst those for which we have ambiguity codes, namely:', \
        ' '.join(ReverseIUPACdict2.keys()) + '. Quitting.', file=sys.stderr)
        raise
  if coverage < MinCovForUpper - 0.5:
    return BaseHere.lower()
  else:
    return BaseHere.upper()

ExpectedBasesNoN = ['A', 'C', 'G', 'T', '-']
NumExpectedBases = len(ExpectedBasesNoN)
def CallEnoughBases(BaseCounts, MinCoverage, coverage):
  '''Analyse base counts to see how many bases need to be called to hit the
  minimum required count.'''

  # Sort the counts from largest to smallest, and sort the associated bases into
  # a matching order.
  SortedBaseCounts, SortedExpectedBases = \
  zip(*sorted(zip(BaseCounts, ExpectedBasesNoN), reverse=True))

  # Iterate through the counts, from largest to smallest, updating the total
  # so far. We should stop once we reach the desired total, but not if the next
  # count is the same as the current one - then we should take that one too.
  # If we reach the end of the list, there's no 'next' to check: we need all the
  # bases.
  CountSoFar = 0  
  for i, count in enumerate(SortedBaseCounts):
    if i == NumExpectedBases:
      NumBasesNeeded = i+1
      break
    CountSoFar += count
    if count == SortedBaseCounts[i+1]:
      continue
    if CountSoFar >= MinCoverage:
      NumBasesNeeded = i+1
      break
  BasesNeeded = SortedExpectedBases[:NumBasesNeeded]
  return CallAmbigBaseIfNeeded(BasesNeeded, coverage)

# Read in the base frequency file
consensus = ''
RefSeq = ''
with open(BaseFreqFile, 'r') as f:

  # Loop through all lines in the file
  for LineNumMin1, line in enumerate(f):

    if LineNumMin1 == 0:
      continue

    # Split up the line into fields separated by commas
    fields = line.split(',')
    if len(fields) != 8:
      print('Line', str(LineNumMin1+1), ',\n', line, 'in', BaseFreqFile, \
      'contains', len(fields), 'fields; expected 8. Quitting', file=sys.stderr)
      exit(1)
    RefBase = fields[1]
    counts = fields[2:]

    # Convert to ints    
    try:
      counts = map(int, counts)
    except ValueError:
      print('Could not understand the base counts as ints on line', \
      str(LineNumMin1+1), ',\n', line, 'in', BaseFreqFile + \
      '. Quitting', file=sys.stderr)
      exit(1)

    # Check positive
    if any(count < 0 for count in counts):
      print('Negative count on line', str(LineNumMin1+1), ',\n', \
      line, 'in', BaseFreqFile + '. Quitting', file=sys.stderr)
      exit(1)

    # Append the reference base
    if len(RefBase) == 1:
      RefSeq += RefBase
    else:
      print('The reference base on line', str(LineNumMin1+1), ',\n', line, \
      'in', BaseFreqFile, 'is', RefBase + \
      '. One character only was expected. Quitting.', file=sys.stderr)
      exit(1)

    # Ignore the count for 'N'. Sum the counts.
    counts = counts[:-1]
    coverage = sum(counts)

    # Call the appropriate character if coverage is below the threshold.
    if coverage < MinCoverage:
      consensus += '?'
      continue

    # Find the base with the highest count (or bases, if they have joint-highest
    # counts).
    MaxCount = max(counts)
    BasesWithMaxCount = [ExpectedBasesNoN[i] for i,count in enumerate(counts) \
    if count == MaxCount]

    # If we're calling the most common base regardless of count:
    if CallMostCommon: 
      BaseToCall = CallAmbigBaseIfNeeded(BasesWithMaxCount, coverage)

    else:
      CountToCallBase = coverage * args.MinFracToCall
      # This next 'if' would also be covered by the 'else', but the explicit
      # 'if' scope is faster and is what is usually needed.
      if MaxCount * len(BasesWithMaxCount) >= CountToCallBase:
        BaseToCall = CallAmbigBaseIfNeeded(BasesWithMaxCount, coverage)
      else:
        BaseToCall = CallEnoughBases(counts, CountToCallBase, coverage)

    consensus += BaseToCall
      
# Replaces gaps that border "no coverage" by "no coverage".
consensus = PropagateNoCoverageChar(consensus)

# Skip positions at which the ref has a gap and the consensus has a gap or
# missing cov.
NewConsensus = ''
NewRefSeq = ''
for ConsensusBase, RefBase in itertools.izip(consensus, RefSeq):
  if RefBase == GapChar and (ConsensusBase == '?' or ConsensusBase == GapChar):
    continue
  NewConsensus += ConsensusBase
  NewRefSeq += RefBase
consensus = NewConsensus
RefSeq = NewRefSeq

ConsensusSeqObj = SeqIO.SeqRecord(Seq.Seq(consensus), \
id=args.consensus_seq_name, description='')
RefSeqObj = SeqIO.SeqRecord(Seq.Seq(RefSeq), \
id=args.ref_seq_name, description='')

SeqIO.write([ConsensusSeqObj, RefSeqObj], sys.stdout, "fasta")

## Thanks Stackoverflow:
#def insert_newlines(string, every=50):
#    lines = []
#    for i in xrange(0, len(string), every):
#        lines.append(string[i:i+every])
#    return '\n'.join(lines)
#print('>' + args.consensus_seq_name)
#print(insert_newlines(consensus))
#print('>' + args.ref_seq_name)
#print(insert_newlines(RefSeq))
