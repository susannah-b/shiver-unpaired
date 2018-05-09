def CalculateReadIdentity(PysamRead, ReferenceSeq):
  '''Calculate the fractional agreement between a read and the ref sequence'''

  positions = PysamRead.get_reference_positions(full_length=True)
  seq = PysamRead.query_sequence.upper()
  ReferenceSeq = ReferenceSeq.upper()
  NumAgreeingBases = 0
  NumDeletions = 0
  LastRefPos = None
  for i, pos in enumerate(positions):
    if pos != None:
      if ReferenceSeq[pos] == seq[i]:
        NumAgreeingBases += 1
      if LastRefPos != None and pos != LastRefPos + 1:
        DeletionSize = pos - LastRefPos - 1
        assert DeletionSize > 0
        NumDeletions += DeletionSize
      LastRefPos = pos
  return float(NumAgreeingBases) / (len(positions) + NumDeletions)

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
