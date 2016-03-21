#!/usr/bin/env python2
from __future__ import print_function

## Author: Chris Wymant, c.wymant@imperial.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''This script removes those columns/positions in alignment
that consist solely of the gap character "-". Output is printed to stdout.'''

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
parser.add_argument('FastaFile', type=File)
args = parser.parse_args()

alignment = AlignIO.read(args.FastaFile, "fasta")
AlignmentLength = alignment.get_alignment_length()
for column in reversed(xrange(AlignmentLength)):
  PureGap = True
  for base in alignment[:, column]:
    if base != '-':
      PureGap = False
      break
  if PureGap:
    alignment = alignment[:, :column] + alignment[:, column+1:]

AlignIO.write(alignment, sys.stdout, 'fasta')
