"""Sequence module defines basic sequence classes for slicable sequences

"""
import sys

_NT2NUM = {'A':2,'a':2,'C':1,'c':1,'G':3,'g':3,'T':0,'t':0}
_NUM2NT = {2:'A',1:'C',3:'G',0:'T'};

class SequenceGeneric:
  """A sequence type that should get overridden
  """
  def __init__(self):
    return
  # must be sliceable
  def __getitem__(self,key):
    sys.stderr.write("must be overriden\n")
    sys.exit()
  def __setitem__(self,key):
    sys.stderr.write("must be overriden\n")
    sys.exit()
  def __delitem__(self,key):
    sys.stderr.write("must be overriden\n")
    sys.exit()
  def __str__(self):
    sys.stderr.write("must be overriden\n")
    sys.exit()
  
class Sequence(SequenceGeneric):
  """The Sequence class in its basic form will run on a string and 
     be very string-like

  :param sequence: take as an input, the sequence
  :type sequence: String
  """
  def __init__(self,sequence):
    self._data = sequence
  def __getitem__(self,key):
    if isinstance(key,slice):
      return Sequence(self._data[key.start:key.stop:key.step])
  def __setitem__(self,key):
    return
  def __delitem__(self,key):
    return
  def __str__(self):
    return self._data
