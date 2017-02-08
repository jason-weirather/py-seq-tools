"""Structure is a module set on defining genomic features

"""
import sys, random, string, uuid
from collections import namedtuple
from seqtools.range import GenomicRange
from seqtools.range.multi import ranges_to_coverage, merge_ranges
from seqtools.sequence import rc
import seqtools.graph
from math import sqrt

MappingGenericOptions = namedtuple('MappingGenericOptions',
   [
    'ref',
    'sequence',
    'name',
    'payload'])
class MappingGeneric(object):
  """Class to describe the mapping of a transcript

     Mapping is sliceable by genomic coordinate. Example:

     All constructors of mapping type have two inputs

     1. rngs 
     2. options

     :param rngs:
     :param options:
     :param options.ref:
     :param options.sequence:
     :param options.name:
     :type rngs: GeneomicRange []
     :type options: namedtuple
     :type options.ref: dict()
     :type options.sequence: String
     :type options.name: String
  """
  def __init__(self,rngs,options=None):
    if not options: options = MappingGeneric.Options()
    self._rngs = rngs
    self._options = options
    self._id = str(uuid.uuid4())

  @staticmethod
  def Options(**kwargs):
      """ A method for declaring options for the class"""
      construct = MappingGenericOptions #IMPORTANT!  Set this
      names = construct._fields
      d = {}
      for name in names: d[name] = None #default values
      """set defaults here"""
      for k,v in kwargs.iteritems():
         if k in names: d[k] = v
         else: raise ValueError('Error '+k+' is not an options property')
      """Create a set of options based on the inputs"""
      return construct(**d)

  @property
  def exons(self):
    """Maybe the most core override this later"""
    return self._rngs

  @property
  def length(self):
    """The mapped length"""
    return sum([x.length for x in self._rngs])

  @property
  def junctions(self):
    """Can be inferred from the exons, this is not implemented yet"""
    sys.stderr.write("Error not implemented junctions\n")
    sys.exit()
    return self._junctions

  def copy(self):
    """A copy of the MappingGeneric

    :return: MappingGeneric copy
    :rtype: MappingGeneric
    """
    return self.__init__(self._rngs,self._options)

  def set_payload(self,val):
    """Set a payload for this object

    :param val: payload to be stored
    :type val: Anything that can be put in a list
    """
    self._options = self._options._replace(payload = val)

  @property
  def payload(self):
    """Get the payload currently being stored

    :return: payload
    :rtype: anything that can be stored in a list
    """
    return self._options.payload

  @property
  def id(self):
    """Return a unique id created for this transcript when it was made

    :return: uuid4 id as a string
    :rtype: string
    """
    return self._id

  def avg_mutual_coverage(self,gpd):
    """get the coverage fraction of each transcript then return the geometric mean

    :param gpd: Another transcript
    :type gpd: Transcript
    :return: avg_coverage
    :rtype: float
    """
    ov = self.overlap_size(gpd)
    if ov == 0: return 0
    xfrac = float(ov) / float(self.get_length())
    yfrac = float(ov) / float(gpd.get_length())
    return sqrt(xfrac*yfrac)

  def overlap_size(self,tx2):
    """Return the number of overlapping base pairs between two transcripts

    :param tx2: Another transcript
    :type tx2: Transcript
    :return: overlap size in base pairs
    :rtype: int
    """
    total = 0
    for e1 in self.exons:
      for e2 in tx2.exons:
        total += e1.overlap_size(e2)
    return total

  def get_exon_count(self):
    """Count the exons in the transcript

    :return: exon count
    :rtype: int
    """
    return len(self.exons)

  def union(self,tx2): # keep direction and name of self
    """Find the union, or perhaps intersection is a better word for it, for two transcripts.  This makes a new transcript.

    :param tx2: transcript 2
    :type tx2: Transcript
    :return: overlapping portion of the transcripts
    :rtype: Transcript
    """
    sys.stderr.write("This is not named or implemented correct as union\n")
    sys.exit()
    self._initialize()
    all = []
    for rng1 in [x.rng for x in self.exons]:
      for rng2 in [y.rng for y in tx2.exons]:
        u = rng1.union(rng2)
        if u: all.append(u)
    if len(all) == 0: return None
    rngs = merge_ranges(all)
    tx = Transcript()
    tx.set_exons_and_junctions_from_ranges(rngs)
    tx._direction = self._direction
    tx._transcript_name = self._transcript_name
    tx._gene_name = self._gene_name
    sys.stderr.write("not updated\n")
    sys.exit()
    return tx

  #def __len__(self):
  #  """Return the length of the mapping in bp. Its the sum of the ranges
  #
  #  :return: length
  #  :rtype: int
  #  """
  #  return sum([x.length for x in self.exons])

  @property
  def sequence(self):
    """A strcutre is defined so get,
    if the sequence is not already there, get the sequence from the reference

    Always is returned on the positive strand for the MappingGeneric

    :param ref_dict: reference dictionary (only necessary if sequence has not been set already)
    :type ref_dict: dict()
    """
    if self._sequence: return self._sequence
    if not self._ref:
      sys.stderr.write("ERROR: sequence is not defined and reference is undefined\n")
      sys.exit()
    self.set_sequence(ref=self._ref)
    return self._sequence

  def set_sequence(self,ref=None,seq=None):
    """use the reference dictionary to set the transcript's sequence

    :param ref: reference dictionary
    :type ref: dict()
    :param seq: sequence
    :type seq: String
    """
    sys.stderr.write("no refactored yet\n")
    sys.exit()
    strand = '+'
    if not self._direction:
      sys.stderr.write("WARNING: no strand information for the transcript\n")
    if self._direction: strand = self._direction
    chr = self.get_chrom()
    seq = ''
    for e in [x.get_range() for x in self.exons]:
      seq += ref_dict[chr][e.start-1:e.end]
    if strand == '-':  seq = rc(seq)
    self._sequence = seq.upper()

  def get_junctions_string(self):
    """Get a string representation of the junctions.  This is almost identical to a previous function.

    :return: string representation of junction
    :rtype: string
    """
    self._initialize()
    return ';'.join([x.get_range_string() for x in self.junctions])

  def junction_overlap(self,tx,tolerance=0):
    """Calculate the junction overlap between two transcripts

    :param tx: Other transcript
    :type tx: Transcript
    :param tolerance: how close to consider two junctions as overlapped (default=0)
    :type tolerance: int
    :return: Junction Overlap Report
    :rtype: Transcript.JunctionOverlap
    """
    self._initialize()
    return JunctionOverlap(self,tx,tolerance)

  def smooth_gaps(self,min_intron):
    """any gaps smaller than min_intron are joined, andreturns a new mapping with gaps smoothed

    :param min_intron: the smallest an intron can be, smaller gaps will be sealed
    :type min_intron: int
    :return: a mapping with small gaps closed
    :rtype: MappingGeneric
    """
    rngs = [self._rngs[0].copy()]
    for i in range(len(self._rngs)-1):
      dist = -1
      if self._rngs[i+1].chr == rngs[-1].chr:
        dist = self._rngs[i+1].start - rngs[-1].end-1
      if dist >= min_intron or dist < 0:
        rngs.append(self._rngs[i+1].copy())
      else:
        rngs[-1].end = self._rngs[i+1].end
    return type(self)(rngs,self._options)
