"""Classes to work with the psl format"""
import seqtools.align
from collections import namedtuple
from seqtools.range import GenomicRange
from seqtools.sequence import rc

PSLFields = namedtuple('PSLFields',
['matches',
'misMatches',
'repMatches',
'nCount',
'qNumInsert',
'qBaseInsert',
'tNumInsert',
'tBaseInsert',
'strand',
'qName',
'qSize',
'qStart',
'qEnd',
'tName',
'tSize',
'tStart',
'tEnd',
'blockCount',
'blockSizes',
'qStarts',
'tStarts'])

PSLOptions = namedtuple('AlignmentOptions',
   ['reference',
    'query_sequence',
    'query_quality'
   ])
class PSL(seqtools.align.Alignment):
  """Class to define a psl line

  :param psl_line: is a psl formated line
  :param reference: a dict/slice accessable sequences
  :param query_sequences: a dict/slice accessable sequences
  :param query_sequence:  just the string that is the query sequence
  :type psl_line: string
  :type refernece: dict()
  :type query_sequences: dict()
  :type query_sequence: string
  """
  def __init__(self,psl_line,options=None):
    if not options: options = PSL.Options()
    self._options = options
    self._line = psl_line.rstrip()
    #self._query_sequence = self._options.query_sequence
    #self._query_quality = self._options.query_quality
    #self._reference = self._options.reference
    self._entries = self._parse_psl_line()
    self._alignment_ranges = None
    self._set_alignment_ranges()

  @property
  def entries(self):
     return self._entries

  @staticmethod
  def Options(**kwargs):
     """ A method for declaring options for the class"""
     construct = PSLOptions #IMPORTANT!  Set this
     names = construct._fields
     d = {}
     for name in names: d[name] = None #default values
     for k,v in kwargs.iteritems():
       if k in names: d[k] = v
       else: raise ValueError('Error '+k+' is not a property of these options')
     """Create a set of options based on the inputs"""
     return construct(**d)

  def __str__(self):
    return self._line

  def get_line(self):
    return self._line

  @property
  def query_sequence(self):
    """Do our overrides parent to get query sequence

    :return: query sequence
    :rtype: string
    """
    return self._options.query_sequence

  @property
  def query_quality(self):
    return self._options.query_quality

  @property
  def reference(self):
    """overrides parent to get the reference genome dict()"""
    return self._options.reference

  #def get_query_length(self):
  #  """overrides parent to get the query length"""
  #  return self.value('qSize')

  def get_PSL(self):
    """Overrides parent to make the PSL generation just return self"""
    return self
  def get_strand(self):
    """same as direction

    :return: strand + or -
    :rtype: char
    """
    return self.value('strand')

  def _parse_psl_line(self):
    f = self._line.rstrip().split("\t")
    if len(f) != 21:
      sys.stderr.write("ERROR: PSL line must contain 21 entries\n")
    return PSLFields(matches=int(f[0]),
       misMatches=int(f[1]),
       repMatches=int(f[2]),
       nCount=int(f[3]),
       qNumInsert=int(f[4]),
       qBaseInsert=int(f[5]),
       tNumInsert=int(f[6]),
       tBaseInsert=f[7],
       strand=f[8], 
       qName=f[9],
       qSize=int(f[10]),
       qStart=int(f[11]),
       qEnd=int(f[12]),
       tName=f[13],
       tSize=int(f[14]),
       tStart=int(f[15]),
       tEnd=int(f[16]),
       blockCount=int(f[17]),
       blockSizes=[int(x) for x in f[18].rstrip(',').split(',')],
       qStarts=[int(x) for x in f[19].rstrip(',').split(',')],
       tStarts=[int(x) for x in f[20].rstrip(',').split(',')])

  # Now we can set things
  # Set list of [target range, query range]
  def _set_alignment_ranges(self):
    self._target_range = GenomicRange(self.entries.tName,self.entries.tStart,self.entries.tEnd)
    self._alignment_ranges = []
    for i in range(0,len(self.entries.blockSizes)):
      trng = GenomicRange(self.entries.tName,self.entries.tStarts[i]+1,self.entries.tStarts[i]+self.entries.blockSizes[i])
      qrng = GenomicRange(self.entries.qName,self.entries.qStarts[i]+1,self.entries.qStarts[i]+self.entries.blockSizes[i])
      self._alignment_ranges.append([trng,qrng])
    return

