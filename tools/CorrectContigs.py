#!/usr/bin/env python
from __future__ import print_function

## Author: Chris Wymant, c.wymant@imperial.ac.uk and Francois Blanquart who
## wrote the original version of this script in R.
## Acknowledgement: we wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''Checks whether contigs appear to need correction by
analysing a comma-separated blast file with variables
qseqid, sseqid, evalue, pident, qlen, qstart, qend, sstart, send
resulting from blasting contigs to a set of reference sequences. By default, we
exit with an error if correction is needed - if contigs have multiple hits
(ignoring any hit fully inside another) or hits are in the reverse direction.
Optionally, the contigs themselves can be supplied and correction is attempted:
cutting the contigs where they have multiple hits and/or reverse-complementing
them where hits are in the reverse direction.
'''

import argparse
import os
import sys
import copy
from Bio import SeqIO
import collections

# Define a function to check files exist, as a type for the argparse.
def File(MyFile):
  if not os.path.isfile(MyFile):
    raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
  return MyFile

# Set up the arguments for this script
parser = argparse.ArgumentParser(description=ExplanatoryMessage)
parser.add_argument('BlastFile', type=File)
parser.add_argument('-C', '--contigs', type=File, \
help='The fasta file of contigs.')
parser.add_argument('-O', '--out-file', help='The file to which corrected' + \
' contigs will be written. If correction is not required, this file will ' +\
'not be created.')
parser.add_argument('-F', '--min-hit-frac', type=float, help='Used to specify'+\
" a minimum fraction of a contig's length that its hit must cover: below this"+\
" we assume correction is needed (in the form of removing the contig or "+\
"deleting the non-blasting segment). This option is ignored if the contigs "+\
"are supplied i.e. if correction is attempted.")
parser.add_argument('-D', '--duplicate-between-hits', action='store_true', \
help='When a contig has multiple hits (after removing any hit fully inside ' +\
'another), we cut the contig. The default behaviour is to cut at the '+\
'midpoint between the two hits (or the middle of the overlap if the hits ' +\
'overlap). With this flag, we duplicate the sequence in between the two hits '+\
'(or the overlap if the hits overlap).')
args = parser.parse_args()
HaveHitFrac = args.min_hit_frac != None

# The --contigs and --out-file flags need each other.
MakeCorrections = args.contigs != None
if MakeCorrections and args.out_file == None:
  print('The --contigs and --out-file flags need each other: use both or', \
  'neither. Quitting.', file=sys.stderr)
  exit(1)

# Check the output file doesn't exist.
if MakeCorrections and os.path.isfile(args.out_file):
  print(args.out_file, 'exists already. Move, rename or delete it, and try', \
  'again. Quitting.', file=sys.stderr)
  exit(1)

# Check the min hit frac is in (0,1)
if HaveHitFrac and not (0 < args.min_hit_frac < 1):
  print('The --min-hit-frac value should be greater than 0 and less than 1.', \
  'Quitting.', file=sys.stderr)
  exit(1)

# columns of the blast output are:
# qseqid means Query Seq-id
# sseqid means Subject Seq-id
# evalue means Expect value
# pident means Percentage of identical matches
# qlen means Query sequence length
# qstart means Start of alignment in query
# qend means End of alignment in query
# sstart means Start of alignment in subject
# send means End of alignment in subject

# Read in the blast hits.
HitDict = collections.OrderedDict()
with open(args.BlastFile) as f:
  for line in f:
    try:
      qseqid, sseqid, evalue, pident, qlen, qstart, qend, sstart, send = \
      line.split(',')
    except ValueError:
      print('The following line in ', args.BlastFile, ' does not contain 9 ', \
      'columns:\n', line, 'Quitting.', sep='', file=sys.stderr)
      exit(1)
    try:
      evalue, pident, qlen, qstart, qend, sstart, send = float(evalue), \
      float(pident), int(qlen), int(qstart), int(qend), int(sstart), int(send)
    except ValueError:
      print('Could not understand columns 3-9 (evalue, pident, qlen, qstart,' ,\
      'qend, sstart, send) on line\n', line, 'in ', args.BlastFile, ' as',\
      'floats & ints. Quitting.', sep='', file=sys.stderr)
      exit(1)
    if qstart >= qend:
      print('qstart greater than or equal to qend on line\n', line, 'in ', \
      args.BlastFile + '. Unexpected. Quitting.', sep='', file=sys.stderr)
      exit(1)
    if sstart == send:
      print('sstart = send on line\n', line, 'in ', args.BlastFile + \
      '. Unexpected. Quitting.', sep='', file=sys.stderr)
      exit(1)
    if min(qstart, sstart, send) < 1:
      print('Non-positive value of qstart, sstart or send on line\n', line, \
      'in ', args.BlastFile + '. Unexpected. Quitting.', sep='', \
      file=sys.stderr)
      exit(1)
    if qend > qlen:
      print('qend greater than qlen on line\n', line, 'in ', args.BlastFile + \
      '. Unexpected. Quitting.', sep='', file=sys.stderr)
      exit(1)
    if qseqid in HitDict:
      HitDict[qseqid].append(\
      [qseqid, sseqid, evalue, pident, qlen, qstart, qend, sstart, send])
    else:
      HitDict[qseqid] = [\
      [qseqid, sseqid, evalue, pident, qlen, qstart, qend, sstart, send]]

# Quit if no hits.
if len(HitDict) == 0:
  print(args.BlastFile, 'contains no hits. Quitting.', file=sys.stderr)
  exit(1)

for contig, hits in HitDict.items():

  # Where a contig has one hit contained entirely inside another, remove the
  # sub-hit.
  if len(hits) > 1:
    SubHitIndices = []
    for i, hit1 in enumerate(hits):
      if i in SubHitIndices:
        continue
      qstart1, qend1 = hit1[5:7]
      qmin1, qmax1 = min(qstart1, qend1), max(qstart1, qend1)
      for j, hit2 in enumerate(hits[i+1:]):
        j = i + j + 1 # offset & zero-based indexing
        if j in SubHitIndices or i in SubHitIndices:
          continue
        qstart2, qend2 = hit2[5:7]
        qmin2, qmax2 = min(qstart2, qend2), max(qstart2, qend2)
        # If the hits start and end at the same place, only one should go:
        if qmin1 == qmin2 and qmax1 == qmax2:
          SubHitIndices.append(j)
        elif qmin1 <= qmin2 and qmax1 >= qmax2:
          SubHitIndices.append(j)
        elif qmin1 >= qmin2 and qmax1 <= qmax2:
          SubHitIndices.append(i)
    assert len(SubHitIndices) == len(set(SubHitIndices)), \
    'Internal error removing sub-hits. Please report to the code author.'
    for i in sorted(SubHitIndices, reverse=True):
      del HitDict[contig][i]
    hits = HitDict[contig]

    # If we're only checking (not correcting) and there are multiple hits or
    # a reverse hit or a too small hit, quit.
    if not MakeCorrections and len(hits) > 1:
      print('Contig correction required (contig', contig, 'has multiple hits,',\
      'after removing any that are fully contained within another). Quitting.',\
      file=sys.stderr)
      exit(1)
  FirstHit = hits[0]
  qseqid, sseqid, evalue, pident, qlen, qstart, qend, sstart, send = FirstHit
  HitFrac = float(qend - qstart + 1) / qlen
  if not MakeCorrections and HaveHitFrac and \
  HitFrac < args.min_hit_frac:
    print('Contig correction required (the hit\n', \
    ' '.join(map(str, FirstHit)), '\nfor contig', contig, 'has a hit fraction',\
    HitFrac, "which is below the specified --min-hit-frac of", \
    str(args.min_hit_frac) + "). Quitting.", file=sys.stderr)
    exit(1)
  if not MakeCorrections and sstart > send:
    print('Contig correction required (the hit\n', \
    ' '.join(map(str, FirstHit)), '\nfor contig', contig + \
    "is in the reverse direction). Quitting.", file=sys.stderr)
    exit(1)

# If we're only checking (not correcting), quit successfully.
if not MakeCorrections:
  print(args.BlastFile, 'seems fine. Quitting successfully.', file=sys.stderr)
  exit(0)

# Read in the contigs
ContigDict = collections.OrderedDict()
for seq in SeqIO.parse(open(args.contigs),'fasta'):
  if seq.id in ContigDict:
    print('Encountered sequence', seq.id, 'a second time in', args.contigs+\
    'Sequence names should be unique. Quitting.', file=sys.stderr)
    exit(1)
  ContigDict[seq.id] = seq

# Check we have a sequence for each hit
UnknownHits = [hit for hit in HitDict.keys() if not hit in ContigDict.keys()]
if len(UnknownHits) != 0:
  print('The following hits in', args.BlastFile, 'do not have a corresponding',\
  'sequence in', args.contigs +':\n', ' '.join(UnknownHits) + \
  '\nQuitting.', file=sys.stderr)
  exit(1)

OutSeqs = []
for ContigName, hits in HitDict.items():

  seq = ContigDict[ContigName]
  SeqLength = len(seq.seq)

  for hit in hits:
    qlen = hit[4]
    if qlen != SeqLength:
      print(args.BlastFile, 'has qlen =', qlen, 'for contig', ContigName, \
      'but that contig has length', SeqLength, 'in', args.contigs + \
      '. Quitting.', file=sys.stderr)
      exit(1)

  NumHits = len(hits)
  if NumHits == 1:
    sstart, send = hits[0][7:9]
    if sstart > send:
      seq.seq = seq.seq.reverse_complement()
    OutSeqs.append(seq)

  else:

    # Sort hits by their start point.
    StartPoints = [hit[5] for hit in hits]
    assert len(StartPoints) == len(set(StartPoints)), \
    'Internal error cutting contigs. Please report to the code author.'
    hits = sorted(hits, key=lambda x:x[3])

    for i, hit in enumerate(hits):

      ThisStart, ThisEnd = hit[5:7]
      if args.duplicate_between_hits:
        if i == 0:
          NextStart, NextEnd = hits[i+1][5:7]
          CutStart = 1
          CutEnd = max(ThisEnd, NextStart)
        elif i == NumHits-1:
          LastStart, LastEnd = hits[i-1][5:7]
          CutStart = min(LastEnd, ThisStart) + 1
          CutEnd = SeqLength
        else:
          LastStart, LastEnd = hits[i-1][5:7]
          NextStart, NextEnd = hits[i+1][5:7]
          CutStart = min(LastEnd, ThisStart) + 1
          CutEnd = max(ThisEnd, NextStart)
      else:
        if i == 0:
          NextStart, NextEnd = hits[i+1][5:7]
          CutStart = 1
          CutEnd = (ThisEnd + NextStart) / 2
        elif i == NumHits-1:
          LastStart, LastEnd = hits[i-1][5:7]
          CutStart = (LastEnd + ThisStart) / 2 + 1
          CutEnd = SeqLength
        else:
          LastStart, LastEnd = hits[i-1][5:7]
          NextStart, NextEnd = hits[i+1][5:7]
          CutStart = (LastEnd + ThisStart) / 2 + 1
          CutEnd = (ThisEnd + NextStart) / 2

      ThisCutSeq = copy.deepcopy(seq)
      ThisCutSeq.id += '.' + str(i+1)
      ThisCutSeq.seq = ThisCutSeq.seq[CutStart-1 : CutEnd]
      
      # Do we need to reverse complement?
      sstart, send = hit[7:9]
      if sstart > send:
        ThisCutSeq.seq = ThisCutSeq.seq.reverse_complement()

      OutSeqs.append(ThisCutSeq)

SeqIO.write(OutSeqs, args.out_file, "fasta")

