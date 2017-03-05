"""Structure is a module set on defining genomic features

"""
import sys, random, string, uuid
from collections import namedtuple
from seqtools.range import GenomicRange, Bed
from seqtools.range.multi import ranges_to_coverage, merge_ranges
from seqtools.sequence import rc
from seqtools.sequence import Sequence
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
    if len(rngs) > 0:
       self._rngs = rngs
    self._options = options
    self._id = str(uuid.uuid4())

  @property
  def options(self): return self._options

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

  def set_reference(self,ref):
     """Assign the reference"""
     self._options = self._options._replace(ref=ref)

  @property
  def reference(self): return self._options.ref

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

  def slice_target(self,chr,start,end):
     """Slice the mapping by the target coordinate
        
        First coordinate is 0-indexed start
        Second coordinate is 1-indexed finish

     """
     # create a range that we are going to intersect with
     trng = Bed(chr,start,end)
     nrngs = []
     for r in self._rngs:
        i = r.intersect(trng)
        if not i: continue
        nrngs.append(i)
     if len(nrngs) == 0: return None
     return MappingGeneric(nrngs,self._options)

  def slice_sequence(self,start,end):
     """Slice the mapping by the position in the sequence
        
        First coordinate is 0-indexed start
        Second coordinate is 1-indexed finish

     """
     #find the sequence length
     l = self.length
     indexstart = start
     indexend = end
     ns = []
     tot = 0
     for r in self._rngs:
        tot += r.length
        n = r.copy()
        if indexstart > r.length:  
           indexstart-=r.length
           continue
        n.start = n.start+indexstart
        if tot > end: 
           diff = tot-end
           n.end -= diff
           tot = end
        indexstart = 0
        ns.append(n)
        if tot == end: break
     if len(ns)==0: return None
     return MappingGeneric(ns,self._options)

  @property
  def sequence(self):
    """A strcutre is defined so get,
    if the sequence is not already there, get the sequence from the reference

    Always is returned on the positive strand for the MappingGeneric

    :param ref_dict: reference dictionary (only necessary if sequence has not been set already)
    :type ref_dict: dict()
    """
    if not self._options.ref:
      raise ValueError("ERROR: sequence is not defined and reference is undefined")
    strand = '+'
    if not self._options.direction:
      sys.stderr.write("WARNING: no strand information for the transcript\n")
    if self._options.direction: strand = self._options.direction
    #chr = self.range.chr
    seq = ''
    for e in [x.range for x in self.exons]:
      seq += str(self._options.ref[e.chr][e.start-1:e.end])
    if strand == '-':  seq = rc(seq)
    return Sequence(seq.upper())

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
