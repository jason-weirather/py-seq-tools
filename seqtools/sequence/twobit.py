"""Sequence module defines basic sequence classes for slicable sequences

"""
import sys
import seqtools.sequence

_NT2NUM = {'A':2,'a':2,'C':1,'c':1,'G':3,'g':3,'T':0,'t':0}
_NUM2NT = {2:'A',1:'C',3:'G',0:'T'};


class Sequence2Bit(seqtools.sequence.SequenceGeneric):
  """For low memory sequence storage
     we are limited to ATCGN for bases
  """
  def __init__(self,sequence=None,base_data=None,mask_data=None,sequence_length=None):
    if base_data and mask_data and sequence_length:
      self._slen = sequence_length
      self._mask = mask_data
      self._data = base_data
    else:
      self._slen = len(str(sequence))
      self._initialize_data()
      self._set_seq(str(sequence)) # set the sequence
  def _set_seq(self,sequence,start=0):
    """Set the sequence from a start position for the length of the sequence"""
    if start+len(sequence) > self._slen: 
      sys.stderr.write("Error not long enough to add\n")
      sys.exit()
    z = 0
    for i in xrange(start, start+len(sequence)):
      self._set_nt(sequence[z],i)
      z+=1

  def _set_nt(self,nuc,pos):
    j = int(pos/4)
    r = pos % 4
    if nuc in _NT2NUM:
      self._data[j] |= _NT2NUM[nuc] << r*2
    # now make the mask
    j = int(pos/8)
    r = pos % 8
    mask = 0
    if nuc not in _NT2NUM: mask = 1
    self._mask[j] |= (mask << r)
  def _get_nt(self,pos):
    j = int(pos/8)
    r = pos % 8
    v = self._mask[j]
    v >>= r
    mask = v & 1
    if mask: return 'N'
    j = int(pos/4)
    r = pos % 4
    v = self._data[j]
    v >>= r*2
    base = _NUM2NT[v & 3];
    return base
  def _initialize_data(self):
    dlen = int(self._slen/4)
    if self._slen % 4 != 0: dlen += 1
    #self._data = c_ubyte*dlen
    self._data = bytearray(dlen)
    nlen = int(self._slen/8)
    if self._slen % 8 != 0: nlen += 1
    #self._mask = c_ubyte*nlen
    self._mask = bytearray(nlen)
    return
  def __str__(self):
    return ''.join([self._get_nt(i) for i in xrange(0, self._slen)])
  def __getitem__(self,key):
    if isinstance(key,slice):
      if key.step:
        sys.stderr.write("Error: step is not useful in a sequence slice\n")
        sys.exit()
      return Sequence2Bit(''.join([self._get_nt(i) for i in xrange(key.start,key.stop)]))
  def __setitem__(self,key):
    return
  def __delitem__(self,key):
    return
