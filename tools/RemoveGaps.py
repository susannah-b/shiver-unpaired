#!/usr/bin/env python2
from __future__ import print_function

## Author: Chris Wymant, c.wymant@imperial.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview: removes gaps from fasta files and justifies them to a consistent
## line length (as required for example by samtools faidx).

import os.path, sys, string

# Check this file is called from the command line with one argument
if len(sys.argv) != 2:
  print('Exactly one argument (a fasta file) should be provided.\nQuitting.', \
  file=sys.stderr)
  exit(1)

FastaFile = sys.argv[1]

# Check that the first argument exists and is a file
if not os.path.isfile(FastaFile):
  print(FastaFile, 'does not exist or is not a file.', file=sys.stderr)
  exit(1)

################################################################################
# USER INPUT
# Any characters that indicate a gap (missing base)
GapChars = '-?'
NumCharsPerLine=50
################################################################################


with open(FastaFile, 'r') as f:

  FoundFirstSequence = False
  AllSequences = []
  for line in f:

    # Strip whitespace & ignore blank lines
    ThisLine = line.strip()
    if ThisLine == '':
      continue

    if ThisLine[0] == '>':
      FoundFirstSequence = True
      NameOfCurrentSequence = ThisLine[1:].strip()
      AllSequences.append([NameOfCurrentSequence,''])
      continue

    if FoundFirstSequence:
      AllSequences[-1][1] += ThisLine

# Check at least one sequence was found
if len(AllSequences) == 0:
  print(FastaFile, 'contains no sequences. This is assumed to be an error.',\
  'Quitting.', file=sys.stderr)
  exit(1)

# Thanks Stackoverflow:
def insert_newlines(string, every=NumCharsPerLine):
    lines = []
    for i in xrange(0, len(string), every):
        lines.append(string[i:i+every])
    return '\n'.join(lines)

for SequenceName,sequence in AllSequences:
  SequenceNoGaps = sequence.translate(string.maketrans("", "", ), GapChars)
  SequenceSplitByNewlines = insert_newlines(SequenceNoGaps)
  print('>'+SequenceName + '\n' + SequenceSplitByNewlines)

