"""Sequence module defines basic sequence classes for slicable sequences

"""
import sys
from collections import namedtuple
from string import maketrans

_NT2NUM = {'A':2,'a':2,'C':1,'c':1,'G':3,'g':3,'T':0,'t':0}
_NUM2NT = {2:'A',1:'C',3:'G',0:'T'};

class SequenceGeneric(object):
  """A sequence type that should get overridden
  """
  def __init__(self,options=None):
    if not options: options = default_options
    self._options = options
    return

  @staticmethod
  def Options(**kwargs):
     """Create a new options namedtuple with only allowed keyword arguments"""
     attributes = ['name','payload']
     Opts = namedtuple('Opts',attributes)
     if not kwargs: return Opts(**dict([(x,None) for x in attributes]))
     kwdict = dict(kwargs)
     for k in [x for x in attributes if x not in kwdict.keys()]: kwdict[k] = None
     return Opts(**kwdict)

  default_options = Options.__func__()

  @property
  def name(self):
    return self._options.name

  @property
  def length(self):
    sys.stderr.write("must be overriden\n")
    sys.exit()
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
  def __init__(self,sequence,options):
    super(Sequence,self).__init__(options)
    self._data = sequence
  def __getitem__(self,key):
    if isinstance(key,slice):
      return Sequence(self._data[key.start:key.stop:key.step],self._options)

  @property
  def length(self):
    return len(self._data)

  @property
  def sequence(self):
    """The nucleotide sequence as a string"""
    return self._data

  def __setitem__(self,key):
    return
  def __delitem__(self,key):
    return
  def __str__(self):
    return self._data

  def rc(self):
    """Reverse complement the sequence and return a Sequence object"""
    return Sequence(rc(self.sequence),self._options)

def rc(seq):
  """Fast reverse complement function using a translation table and slice

  :param seq: string to reverse complement
  :type seq: string
  :return: reverse complemented sequence
  :rtype: string
  """
  complement = maketrans('ACTGUNXactgunx','TGACANXtgacanx')
  return seq.translate(complement)[::-1]
