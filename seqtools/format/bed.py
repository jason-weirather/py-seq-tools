import sys, uuid
from collections import namedtuple
import seqtools.structure.transcript
from seqtools.range import GenomicRange

Bed12Options = namedtuple('Bed12Options',
   ['sequence',
    'ref',
    'gene_name',
    'payload'])

Bed12Fields = namedtuple('Bed12Fields',
   [
    'chrom',
    'chromStart',
    'chromEnd',
    'name',
    'score',
    'strand',
    'thickStart',
    'thickEnd',
    'itemRgb',
    'blockCount',
    'blockSizes',
    'blockStarts'])

class Bed12(seqtools.structure.transcript.Transcript):
  """ Bed format with 9 optional fields 

  :param bed_line: one line of a bed file
  :type bed_line: string
  """
  def __init__(self,bed_line,options=None):
    if not options: options = Bed12.Options()
    self.entries = self._line_to_entry(bed_line)
    rngs = self._get_ranges()
    super(Bed12,self).__init__(rngs,
       super(Bed12,self).Options(
          #fill in some options of Transcript type we are overriding
          direction = self.entries.strand,
          ref = options.ref,
          name = self.entries.name,
          gene_name = options.gene_name,
          payload = options.payload
       ))
    self._line = bed_line.rstrip()

  @property
  def bed12_line(self):
     return self._line

  @staticmethod
  def Options(**kwargs):
      """ A method for declaring options for the class"""
      construct = Bed12Options #IMPORTANT!  Set this
      names = construct._fields
      d = {}
      for name in names: d[name] = None #default values
      """set defaults here"""
      for k,v in kwargs.iteritems():
         if k in names: d[k] = v
         else: raise ValueError('Error '+k+' is not an options property')
      """Create a set of options based on the inputs"""
      return construct(**d)


  def _get_ranges(self):
    rngs = []
    for i in range(0,self.entries.blockCount):
      rng = GenomicRange(self.entries.chrom,\
            self.entries.chromStart+self.entries.blockStarts[i]+1,\
            self.entries.chromStart+self.entries.blockStarts[i]+self.entries.blockSizes[i])
      rngs.append(rng)
    return rngs

  def __str__(self):
    return self.get_bed_line()  

  def _line_to_entry(self,line):
    """parse the line into entries and keys"""
    f = line.rstrip().split("\t")
    """
    'chrom'
    'chromStart'
    'chromEnd'
    'name'
    'score'
    'strand'
    'thickStart'
    'thickEnd'
    'itemRgb'
    'blockCount'
    'blockSizes'
    'blockStarts'
    """
    return Bed12Fields(
       f[0],
       int(f[1]),
       int(f[2]),
       f[3],
       int(f[4]),
       f[5],
       int(f[6]),
       int(f[7]),
       [int(x) for x in f[8].rstrip(',').split(',')],
       int(f[9]),
       [int(x) for x in f[10].rstrip(',').split(',')],
       [int(x) for x in f[11].rstrip(',').split(',')])

