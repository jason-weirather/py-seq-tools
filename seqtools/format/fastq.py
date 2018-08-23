import re, sys
import seqtools.sequence
from collections import namedtuple

class FASTQStream:
  """Iterable Stream"""
  def __init__(self,fh):
    self.fh = fh
  def __iter__(self):
    return self
  def __next__(self):
    return self.next()
  def next(self):
    v = self.get_entry()
    if not v:
      raise StopIteration
    else:
      return v
  def get_entry(self):
    line1 = self.fh.readline().rstrip()
    if not line1: return None
    line2 = self.fh.readline().rstrip()
    if not line2: return None
    line3 = self.fh.readline().rstrip()
    if not line3: return None
    line4 = self.fh.readline().rstrip()
    if not line4: return None
    return FASTQ("\n".join([line1,line2,line3,line4]))

FASTQOptions = namedtuple('FASTQOptions',
   ['name',
    'header',
    'payload'])

class FASTQ(seqtools.sequence.Sequence):
  """ fastq single entry

  :param v: one entry
  :type v: string
  """
  def __init__(self,v):
    self.lines = v.rstrip().split("\n")
    header = self.lines[0][1:]
    name = re.match('(\S+)',header).group(1)
    seq = self.lines[1]
    opts = FASTQ.Options(name=name,header=header)
    super(FASTQ,self).__init__(seq,opts)
    self.qual = self.lines[3]

  @staticmethod
  def Options(**kwargs):
      """ A method for declaring options for the class"""
      construct = FASTQOptions #IMPORTANT!  Set this
      names = construct._fields
      d = {}
      for name in names: d[name] = None #default values
      """set defaults here"""
      for k,v in kwargs.items():
         if k in names: d[k] = v
         else: raise ValueError('Error '+k+' is not an options property')
      """Create a set of options based on the inputs"""
      return construct(**d)

  default_options = Options.__func__()

  def __getitem__(self,key):
    if isinstance(key,slice):
      newseq = str(self.sequence)[key.start:min(key.stop,len(self.sequence))]
      newqual = self.qual[key.start:min(key.stop,len(self.sequence))]
      return FASTQ('@'+self.header+"\n"+newseq+"\n"+self.lines[2]+"\n"+newqual)
    return {'name':self.name,'seq':str(self.sequence),'qual':self.lines[3]}[key]

  @property
  def name(self): return self._options.name
  @property
  def header(self): return self._options.header

  def rc(self):
    return FASTQ('@'+self.header+"\n"+seqtools.sequence.rc(self.sequence)+"\n"+self.lines[2]+"\n"+self.qual[::-1]+"\n")
  def copy(self):
    return FASTQ(self.FASTQ())
  def FASTQ(self):
    return "\n".join(self.lines)+"\n"
  def __str__(self):
    return self.FASTQ()
