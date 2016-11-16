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
