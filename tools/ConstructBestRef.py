#!/usr/bin/env python2
from __future__ import print_function

## Author: Chris Wymant, c.wymant@imperial.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251 
##
## Overview: we construct the optimal reference for a given sample by flattening
## its (de novo assembled) contigs with a set of references.
## How it works:
## * We flatten the contigs, taking the base of the longest one where they
## disagree (expecting that the de novo assembly makes the longest contig first
## with the dominant i.e. most numerous set of short reads).
## * We compare every reference to the flattened contigs and count the number of
## agreeing bases. The comparison is not done in the gaps between contigs nor
## before/after the reference/contigs. When we are inside a contig and inside
## the reference, however, gaps count - here they are taken to represent a
## genuine deletion, so the gap character holds equal footing with a base.
## * Starting with the reference with the highest score, we elongate it both
## directions using references or progressively lower scores. We stop when both
## edges reach the edges of the alignment, or when a reference is reached which
## has a score of zero. This defines the elongated best reference.
## * We fill in gaps before, between and after (but not within) the contigs
## using the elongated best reference. This defines the constructed best
## reference.
##
################################################################################
## USER INPUT
# The character that indicates a gap (missing base). If there are multiple such
# characters, put the one to be used in the output first.
GapChars = '-.?'
# Excise unique insertions present in the contigs but not any of the references?
ExciseUniqueInsertions = False
################################################################################

FloatComparisonTolerance = 1e-6

import os.path, sys, collections
from string import maketrans
from AuxiliaryFunctions import ReadSequencesFromFile

# Check this file is called from the command line with the correct number of
# arguments.
if len(sys.argv) < 3:
  print('Incorrect number of arguments given.\n'+\
  'Usage:\n', sys.argv[0], 'AlignmentFile NameOfContig1',\
  '[NameOfContig2...]', file=sys.stderr)
  exit(1)
AlignmentFile = sys.argv[1]
ContigNames = sys.argv[2:]

# Check input data files exist and are files
if not os.path.isfile(AlignmentFile):
  print(AlignmentFile, 'does not exist or is not a file.', file=sys.stderr)
  exit(1)

# Check all contig names are unique
CounterObject = collections.Counter(ContigNames)
DuplicatedContigNames = [i for i in CounterObject if CounterObject[i]>1]
if len(DuplicatedContigNames) != 0:
  for ContigName in DuplicatedContigNames:
    print('Contig name', ContigName, 'was duplicated in the arguments.', \
    file=sys.stderr)
  print('All contig names should be unique. Exiting.', file=sys.stderr)
  exit(1)

# Read in the sequences from the alignment file (into a dictionary)
AllSeqsDict, AlignmentLength = ReadSequencesFromFile(AlignmentFile)

# Separate sequences into references and contigs
RefDict = {}
ContigDict = {}
for SeqName in AllSeqsDict:
  if SeqName in ContigNames:
    ContigDict[SeqName] = AllSeqsDict[SeqName]
  else:
    RefDict[SeqName] = AllSeqsDict[SeqName]

# Check we found all the contigs
MissingContigs = [contig for contig in ContigNames if not contig in ContigDict]
if len(MissingContigs) != 0:
  for ContigName in MissingContigs:
    print('Contig "'+ContigName+'" was not found in the alignment file.', \
    file=sys.stderr)
  print('Exiting.', file=sys.stderr)
  exit(1)

# Ensure a unique character for gaps
GapChar = GapChars[0]
if len(GapChars) > 1:
  OtherGapChars = GapChars[1:]
  GapCharRepeat = GapChar*len(OtherGapChars)
  for SeqName,seq in ContigDict.items():
    ContigDict[SeqName] = seq.translate(maketrans(OtherGapChars,GapCharRepeat))
  for SeqName,seq in RefDict.items():
    RefDict[SeqName] = seq.translate(maketrans(OtherGapChars,GapCharRepeat))

# A function we'll need more than once:
def FindSeqStartAndEnd(SeqName,seq):
  '''Find the 0-based positions of the start and end of the sequence.'''
  StartOfSeq = 0
  try:
    while seq[StartOfSeq] == GapChar:
      StartOfSeq += 1
  except IndexError:
    print(SeqName, "has no bases - it's just one big gap.\nQuitting.", \
    file=sys.stderr)
    exit(1)
  EndOfSeq = AlignmentLength-1
  while seq[EndOfSeq] == GapChar:
    EndOfSeq -= 1
  return StartOfSeq,EndOfSeq

# Find the start, end, length and gap fraction of each contig.
# Find the start of the left-most ('first') contig and the end of the right-most
# ('last') contig.
ContigStartsAndEnds = {}
ContigLengths = {}
ContigGapFractions = {}
TotalGapsInContigs = 0
for ContigName,ContigSeq in ContigDict.items():
  StartOfContig, EndOfContig = FindSeqStartAndEnd(ContigName, ContigSeq)
  ContigStartsAndEnds[ContigName] = [StartOfContig,EndOfContig]
  NumGapsInAlignedContig = ContigSeq.count(GapChar)
  NumInternalGaps = ContigSeq[StartOfContig:EndOfContig+1].count(GapChar)
  TotalGapsInContigs += NumInternalGaps
  ContigLengths[ContigName] = AlignmentLength -NumGapsInAlignedContig
  ContigGapFractions[ContigName] = \
  float(NumInternalGaps)/(EndOfContig+1 - StartOfContig)

#print(' '.join(map(str, sorted(ContigGapFractions.values(), reverse=True) )))
#print(float(TotalGapsInContigs)/sum(ContigLengths.values()))
#exit(0)

AllContigStarts = []
AllContigEnds = []
for [ContigStart,ContigEnd] in ContigStartsAndEnds.values():
  AllContigStarts.append(ContigStart)
  AllContigEnds.append(ContigEnd)
StartOfFirstContig = min(AllContigStarts)
EndOfLastContig = max(AllContigEnds)

GappiestContigName, GapFraction = \
sorted(ContigGapFractions.items(), key=lambda x:x[1], reverse=True)[0]

# Flatten the contigs.
# Repeat the gap character until the start of the first contig.
# Then at each position: if no contig has a base there, call a gap; if only one
# contig has a base there, use that base; if multiple contigs have bases there,
# use the base of the longest contig. Then repeat the gap character until the
# end.
FlattenedContigsSeq = GapChar * StartOfFirstContig
for position in range(StartOfFirstContig,EndOfLastContig+1):
  DictOfBasesHere = {}
  for ContigName,ContigSeq in ContigDict.items():
    base = ContigSeq[position]
    if base != GapChar:
      DictOfBasesHere[ContigName] = ContigSeq[position]
  if len(DictOfBasesHere) == 0:
    BaseHere = GapChar
  elif len(DictOfBasesHere) == 1:
    BaseHere = DictOfBasesHere.values()[0]
  else:
    LengthOfLongestContigWithBaseHere = 0
    for ContigName in DictOfBasesHere:
      if ContigLengths[ContigName] > LengthOfLongestContigWithBaseHere:
        LongestContigWithBaseHere = ContigName
        LengthOfLongestContigWithBaseHere = ContigLengths[ContigName]
    BaseHere = DictOfBasesHere[LongestContigWithBaseHere]
  FlattenedContigsSeq += BaseHere
FlattenedContigsSeq += GapChar * (AlignmentLength - EndOfLastContig -1)

# Make a list, of the same length of the alignment, of integers: each one
# counting the number of contigs with coverage there. Gaps inside contigs get
# counted as coverage; gaps between contigs get a count of 0.
ContigCoverageByPosition = [0 for n in range(0,AlignmentLength)]
for [start,end] in ContigStartsAndEnds.values():
  for position in range(start,end+1):
    ContigCoverageByPosition[position] += 1

#TotalContigCoverage = sum([1 for base in FlattenedContigsSeq if base != GapChar])
#sys.stdout.write(str(len(ContigDict)) +' '+ str(TotalContigCoverage) +' ')
#exit(0)

# Count the number of positions where exactly one contig has a base and no
# reference has a base - 'unique insertions'. Replace such positions by gaps if
# desired.
# Count how many deletions there are inside contigs.
FlattenedContigsSeq_NoNewInsertions = ''
LengthOfUniqueInsertions = 0
LengthOfDeletions = 0
for position in range(0,AlignmentLength):
  NumContigsHere = ContigCoverageByPosition[position]
  ContigBase = FlattenedContigsSeq[position]
  if NumContigsHere > 0 and ContigBase == GapChar:
    LengthOfDeletions += 1
  if NumContigsHere == 1:
    NoRefHasBaseHere = True
    for RefName,RefSeq in RefDict.items():
      if RefSeq[position] != GapChar:
        NoRefHasBaseHere = False
        break
    if NoRefHasBaseHere:
      LengthOfUniqueInsertions += 1
    if ExciseUniqueInsertions and NoRefHasBaseHere:
      FlattenedContigsSeq_NoNewInsertions += GapChar
    else:
      FlattenedContigsSeq_NoNewInsertions += ContigBase
  else:
    FlattenedContigsSeq_NoNewInsertions += ContigBase
if ExciseUniqueInsertions:
  FlattenedContigsSeq = FlattenedContigsSeq_NoNewInsertions

# For each reference, find the fraction of positions with both reference and
# contig coverage where the two are in agreement. Record the reference start and
# end.
ListOfRefsAndScores = []
for RefName,RefSeq in RefDict.items():
  NumBasesAgreeing = 0
  OverlapLength = 0
  StartOfRef, EndOfRef = FindSeqStartAndEnd(RefName, RefSeq)
  for position in range(StartOfRef, EndOfRef+1):
    if ContigCoverageByPosition[position] > 0:
      OverlapLength += 1
      if RefSeq[position] == FlattenedContigsSeq[position]:
        NumBasesAgreeing += 1
  if OverlapLength == 0:
    FractionalAgreement = 0
  else:
    FractionalAgreement = float(NumBasesAgreeing)/OverlapLength
  ListOfRefsAndScores.append([RefName, StartOfRef, EndOfRef, \
  FractionalAgreement])
  #NumBasesAgreeing])

# Check at least one reference matches at least one position!
if all(item[3] < FloatComparisonTolerance for item in ListOfRefsAndScores):
  print('No reference matches the contigs at any position! We assume this is',\
  'is an error.\nQuitting.', file=sys.stderr)
  exit(1)

# Sort the references - the ones closest to the contigs first.
ListOfRefsAndScores = sorted(ListOfRefsAndScores, key=lambda x:x[3], \
reverse=True)

# Start with the best ref, and iteratively extend it using longer references
# with lower scores. Break if we reach a reference with zero score, or if our
# construct becomes as long as the alignment.
BestRefName, BestRefStart, BestRefEnd, BestRefScore = ListOfRefsAndScores[0]
ElongatedRef = RefDict[BestRefName]
ElongatedRefStart, ElongatedRefEnd = BestRefStart, BestRefEnd
for RefName, StartOfRef, EndOfRef, NumBasesAgreeing in ListOfRefsAndScores[1:]:
  if NumBasesAgreeing == 0 or \
  (ElongatedRefStart == 0 and ElongatedRefEnd == AlignmentLength-1):
    break
  ThisRef = RefDict[RefName]
  if StartOfRef < ElongatedRefStart:
    ElongatedRef = \
    ThisRef[:ElongatedRefStart] + ElongatedRef[ElongatedRefStart:]
    ElongatedRefStart = StartOfRef
  if EndOfRef > ElongatedRefEnd:
    ElongatedRef = \
    ElongatedRef[:ElongatedRefEnd+1] + ThisRef[ElongatedRefEnd+1:]
    ElongatedRefEnd = EndOfRef

# Fill in gaps in contig coverage using the elongated best reference from above.
ConstructedRef = ''
for position in range(0,AlignmentLength):
  if ContigCoverageByPosition[position] > 0:
    ConstructedRef += FlattenedContigsSeq[position]
  else:
    ConstructedRef += ElongatedRef[position]

# Inserting line breaks: thanks Stackoverflow:
def insert_newlines(string, every=50):
  lines = []
  for i in xrange(0, len(string), every):
    lines.append(string[i:i+every])
  return '\n'.join(lines)

# Print output.
#TotalContigCoverage = \
#sum([1 for coverage in ContigCoverageByPosition if coverage > 0])
print('>ContigsFlattenedWith_'+BestRefName)
print(insert_newlines(ConstructedRef))
print('>'+BestRefName+'_elongated')
print(insert_newlines(ElongatedRef))

