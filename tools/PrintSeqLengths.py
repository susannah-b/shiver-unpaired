#!/usr/bin/env python2
from __future__ import print_function

## Author: Chris Wymant, c.wymant@imperial.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''This script takes one fasta file as an argument and
prints the length of each sequence found therein, optionally ignoring gaps.'''

import argparse
import os
import sys
import string
from Bio import SeqIO

# Define a function to check files exist, as a type for the argparse.
def File(MyFile):
  if not os.path.isfile(MyFile):
    raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
  return MyFile

# Set up the arguments for this script
ExplanatoryMessage = ExplanatoryMessage.replace('\n', ' ').replace('  ', ' ')
parser = argparse.ArgumentParser(description=ExplanatoryMessage)
parser.add_argument('FastaFile', type=File)
parser.add_argument('-g', '--include-gaps', action='store_true', \
help='include gap characters, "-" and "?", ignored by default')
parser.add_argument('-C', '--ignore-lower-case', action='store_true', \
help="Doesn't count lower case letters.")
parser.add_argument('-1', '--first-seq-only', action='store_true', \
help='''Print the length of the first sequence only. Faster for alignments
containing many sequences (in which case you'll likely want the -g option
too!).''')
parser.add_argument('-F', '--fragments', action='store_true', help='''Print the
length of the 'fragments' of the sequence, i.e. those parts that are separated
by one or more '?' characters (denoting missing sequence data). If a sequence
consists of N fragments, we will print N values followed by a zero (meaning that
the length of the (N+1)th fragment is zero).''')
args = parser.parse_args()

# Can't use --include-gaps and --fragments
if args.include_gaps and args.fragments:
  print('The --include-gaps and --fragments options cannot both be used at',
  'once. Quitting.', file=sys.stderr)
  exit(1)

SeqLengths = []
for seq in SeqIO.parse(open(args.FastaFile),'fasta'):
  if not args.include_gaps:
    seq.seq = seq.seq.ungap("-")
    if not args.fragments:
      seq.seq = seq.seq.ungap("?")
  if args.ignore_lower_case:
    seq.seq = ''.join(x for x in seq.seq if not x.islower())
  if args.fragments:
    FragLengths = [len(frag) for frag in seq.seq.split('?') if len(frag) > 0]
    FragLengths = sorted(FragLengths, reverse=True)
    FragLengths.append(0)
    SeqLengths.append([seq.id] + FragLengths)
  else:
    SeqLengths.append([seq.id, len(seq.seq)])
  if args.first_seq_only:
    break

for data in SeqLengths:
  print(' '.join(map(str,data)))

