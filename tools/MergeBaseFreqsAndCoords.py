#!/usr/bin/env python
from __future__ import print_function

## Author: Chris Wymant, c.wymant@imperial.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''TODO'''

import sys
import os
import csv
import argparse

# Define a function to check files exist, as a type for the argparse.
def File(MyFile):
  if not os.path.isfile(MyFile):
    raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
  return MyFile

# Set up the arguments for this script
parser = argparse.ArgumentParser(description=ExplanatoryMessage)
parser.add_argument('BaseFreqFile', type=File, help='''A file of base
frequencies, in the format output by shiver (more specifically, by the script
shiver/tools/AnalysePileup.py which is run as part of shiver).''')
parser.add_argument('CoordsFile', type=File, help='''A file detailing a
coordinate tranformation, in the format output by shiver (more specifically, by
the script shiver/tools/MergeAlignments.py which is run as part of shiver).''')
args = parser.parse_args()

RefPosColumnInBaseFreqFile = 1
RefPosColumnInCoordsFile = 2
AlnPosColumnInCoordsFile = 1

RefToAlnPosDict = {}
with open(args.CoordsFile, 'r') as f:
  reader = csv.reader(f, delimiter=',') #, quotechar='"')
  for LineNumberMin1, fields in enumerate(reader):
    if LineNumberMin1 == 0:
      continue
    try:
      RefPos = fields[RefPosColumnInCoordsFile-1]
      AlnPos = fields[AlnPosColumnInCoordsFile-1]
    except IndexError:
      print('Not enough fields on line', LineNumberMin1+1, 'in', \
      args.CoordsFile +'. Quitting.', file=sys.stderr)
      exit(1)
    if RefPos != '-':
      RefToAlnPosDict[RefPos] = AlnPos

with open(args.BaseFreqFile, 'r') as f:
  for LineNumberMin1, line in enumerate(f):
    if LineNumberMin1 == 0:
      OutputString = 'Position in alignment,' + line
      continue
    fields = line.split(',')
    try:
      RefPos = fields[RefPosColumnInBaseFreqFile-1]
    except IndexError:
      print('Not enough fields on line', LineNumberMin1+1, 'in', \
      args.BaseFreqFile +'. Quitting.', file=sys.stderr)
      exit(1)
    try:  
      AlnPos = RefToAlnPosDict[RefPos]
    except KeyError:
      AlnPos = '-'
    OutputString += str(AlnPos) + ',' +line

OutputString = OutputString.rstrip()

print(OutputString)
