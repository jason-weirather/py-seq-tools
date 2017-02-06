"""a module to work with higher level transcript objects

"""
import sys, random, string, uuid
from collections import namedtuple
from seqtools.range import GenomicRange
from seqtools.range.multi import ranges_to_coverage, merge_ranges
from seqtools.sequence import rc
import seqtools.graph
from math import sqrt
from seqtools.structure.transcript.converters import transcript_to_gpd_line, transcript_to_fake_psl_line


class Transcript(seqtools.structure.MappingGeneric):
  """Class to describe a transcript

     This is a basic transcript where all the exons are in order and on the same chromosome

     Mapping is sliceable by genomic coordinate. Example:

     * mymapping - a mapping that ranges between 1 and 100
     * mymapping[1:99] - a new mapping that rangse between 2 and 99
     
     :param rngs:
     :param options:
     :param options.direction:
     :param options.ref:
     :param options.sequence:
     :param options.name:
     :param options.gene_name:
     :type rngs: GeneomicRange []
     :type options: namedtuple
     :type options.direction: Char
     :type options.ref: dict()
     :type options.sequence: String
     :type options.name: String
     :type options.gene_name: String     
  """

  def __init__(self,rngs,options=None):
    if not options: options = {'direction':None,
         'ref':None,
         'sequence':None,
         'name':None,
         'gene_name':None,
         'payload':None}
    super(Transcript,self).__init__(rngs,options)
    self._rngs = rngs
    self._options = options
  

  def __getitem__(self,key):
    """To handle slicing the transcript"""
    if isinstance(key,slice):
      if key.step:
        sys.stderr.write("Step is not apporipriate for slicing a Mapping\n")
      trimmed = trim_ordered_range_list(self._rngs,key.start)
      return self.__init__(trimmed,self._options)

  def __setitem__(self,key):
    return
  def __delitem__(self,key):
    return

  @property
  def range(self):
    """Get the range from the leftmost exon to the rightmost

    :return: total range
    :rtype: GenomicRange
    """
    return GenomicRange(self._rngs[0].chr,self._rngs[0].start,self._rngs[-1].end)

  #def get_range(self):
  #  return self.range

  def set_strand(self,dir):
    """Set the strand (direction)

    :param dir: direction + or -
    :type dir: char
    """
    self._direction = dir

  @property
  def strand(self):
    """Get the strand

    :return: direction + or -
    :rtype: char
    """
    return self._direction

  @property
  def direction(self):
    """alias for strand"""
    return self._direction

  @property
  def chr(self):
    """the reference chromosome. greedy return the first chromosome in exon array

    :return: chromosome
    :rtype: string
    """
    if len(self.exons)==0: 
      sys.stderr.write("WARNING can't return chromsome with nothing here\n")
      return None
    return self._rngs[0].chr

  @property
  def junctions(self):
    """Can be inferred from the exons, this is not implemented yet"""
    sys.stderr.write("Error not implemented junctions\n")
    sys.exit()
    return self._junctions

  def get_gpd_line(self,transcript_name=None,gene_name=None,strand=None):
    """Get the genpred format string representation of the mapping"""
    return transcript_to_gpd_line(self,transcript_name,gene_name,strand)

  def set_gene_name(self,name):
    """assign a gene name

    :param name: name
    :type name: string
    """
    self._options = self._options._replace(gene_name = name)

  def get_gene_name(self):
    """retrieve the gene name

    :return: gene name
    :rtype: string
    """
    self._options.gene_name
    return self._gene_name
  def set_transcript_name(self,name):
    """assign a transcript name

    :param name: name
    :type name: string
    """
    self._options = self._options._replace(name = name)

  def get_transcript_name(self):
    """retrieve the transcript name

    :return: transcript name
    :rtype: string
    """
    return self._options.name

  def get_fake_psl_line(self,ref):
    """Convert a mapping to a fake PSL line"""
    return transcript_to_fake_psl_line(self,ref)

  def get_junctions_string(self):
    """Get a string representation of the junctions.  This is almost identical to a previous function.

    :return: string representation of junction
    :rtype: string
    """
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
    return JunctionOverlap(self,tx,tolerance)

  def exon_overlap(self,tx,multi_minover=10,multi_endfrac=0,multi_midfrac=0.8,single_minover=50,single_frac=0.5,multi_consec=True):
    """Get a report on how mucht the exons overlap

    :param tx:
    :param multi_minover: multi-exons need to overlap by at lest this much to be considered overlapped (default 10)
    :param multi_endfrac: multi-exons need an end fraction coverage of at least this by default (default 0)
    :param multi_midfrac: multi-exons need (default 0.8) mutual coverage for internal exons
    :parma single_minover: single-exons need at least this much shared overlap (default 50)
    :param single_frac: at least this fraction of single exons must overlap (default 0.5)
    :parma multi_consec: exons need to have multiexon consecutive mapping to consider it a match (default True)
    :type tx:
    :type multi_minover: int
    :type multi_endfrac: float
    :type multi_midfrac: float
    :type single_minover: int
    :type single_frac: float
    :type multi_consec: bool
    :return: ExonOverlap report
    :rtype: Transcript.ExonOverlap
    """
    return Transcript.ExonOverlap(self,tx,multi_minover,multi_endfrac,multi_midfrac,single_minover,single_frac,multi_consec=multi_consec)

class ExonOverlap:
    """class to describe exon overlap

    :param tx:
    :param multi_minover: multi-exons need to overlap by at lest this much to be considered overlapped (default 10)
    :param multi_endfrac: multi-exons need an end fraction coverage of at least this by default (default 0)
    :param multi_midfrac: multi-exons need (default 0.8) mutual coverage for internal exons
    :parma single_minover: single-exons need at least this much shared overlap (default 50)
    :param single_frac: at least this fraction of single exons must overlap (default 0.5)
    :parma multi_consec: exons need to have multiexon consecutive mapping to consider it a match (default True)
    :type tx:
    :type multi_minover: int
    :type multi_endfrac: float
    :type multi_midfrac: float
    :type single_minover: int
    :type single_frac: float
    :type multi_consec: bool
    :return: ExonOverlap report
    :rtype: Transcript.ExonOverlap
    """
    def __init__(self1,tx_obj1,tx_obj2,multi_minover=10,multi_endfrac=0,multi_midfrac=0.8,single_minover=50,single_frac=0.5,multi_consec=True):
      self1.tx_obj1 = tx_obj1
      self1.tx_obj2 = tx_obj2
      self1.multi_minover = multi_minover # multi-exon minimum overlap of each exon
      self1.multi_endfrac = multi_endfrac # multi-exon minimum fractional overlap of first or last exon
      self1.multi_midfrac = multi_midfrac # multi-exon minimum fractional overlap of internal exons
      self1.multi_consec = multi_consec # require consecutive exons for exon overlap of multi_exon
      self1.single_minover = single_minover # single-exon minimum overlap
      self1.single_frac = single_frac #single-exon minimum overlap
      self1.overs = [] # set by calculate_overlap()
      #self1.dif1 = []
      #self1.dif2 = []
      self1.calculate_overlap()
      if len(self1.overs) == 0: return None# nothing to analyze
      if self1.tx_obj1.get_exon_count() > 1 and self1.tx_obj1.get_exon_count() > 1 \
         and self1.multi_consec and len(self1.overs) < 2: 
        return None #not enough to consider multi exon transcript overlap
      self1.analyze_overs()
      if self1.tx_obj1.get_exon_count() > 1 and self1.tx_obj1.get_exon_count() > 1 \
        and self1.multi_consec and (min(self1.dif1) != 1 or min(self1.dif2) !=1): 
        return None #not enough to consider multi exon transcript overlap
      
    def __nonzero__(self1):
      if len(self1.overs) > 0: return True
      return False
  
    def match_exon_count(self1):
      """Total number of exons that overlap

      :return: matched exon count
      :rtype: int
      """
      return len(self1.overs)

    def consecutive_exon_count(self1):
      """Best number of consecutive exons that overlap

      :return: matched consecutive exon count
      :rtype: int
      """
      best = 1
      consec = 1
      for i in range(0,len(self1.dif1)):
        if self1.dif1[i] == 1 and self1.dif2[i] == 1:
          consec += 1
        else: 
          consec = 1
        if consec > best:
          best = consec
      return best

    def is_subset(self1):
      """ Return value if tx_obj2 is a complete subset of tx_obj1 or 
          tx_obj1 is a complete subset of tx_obj2

          Values are:

          * Return 1: Full overlap (mutual subests) 
          * Return 2: two is a subset of one
          * Return 3: one is a subset of two
          * Return False if neither is a subset of the other

      :return: subset value
      :rtype: int
      """
      if len(self1.overs) == 0: return False
      if len(self1.dif1) > 0: # make sure they are consecutive if more than one
        if max(self1.dif1) != 1 or max(self1.dif2) != 1: return False
      onecov = self1.start1 and self1.end1
      twocov = self1.start2 and self1.end2
      if onecov and twocov:
        return 1
      elif twocov: return 2
      elif onecov: return 3
      return False

    def is_full_overlap(self1):
      """true if they are a full overlap

      :return: is full overlap
      :rtype: bool
      """
      if len(self1.overs) == 0: return False
      if len(self1.dif1) > 0:
        if max(self1.dif1) != 1 or max(self1.dif2) != 1: return False
      if self1.start1 and self1.end1 and self1.start2 and self1.end2:
        return True
      return False

    def is_compatible(self1):
      """ Return True if the transcripts can be combined together 

      :return: can be combined together
      :rtype: bool
      """
      if len(self1.overs) == 0: return False
      if len(self1.dif1) > 0:
        if max(self1.dif1) != 1 or max(self1.dif2) != 1: return False
      # If we are still here it is a single run
      if (self1.start1 or self1.start2) and (self1.end1 or self1.end2):
        return True
      return False

    def analyze_overs(self1):
      """A helper function that prepares overlap and consecutive matches data"""
      #check for full overlap first
      self1.dif1 = [self1.overs[i][0]-self1.overs[i-1][0] for i in range(1,len(self1.overs))]
      self1.dif2 = [self1.overs[i][1]-self1.overs[i-1][1] for i in range(1,len(self1.overs))]
      #see if it starts and ends on first or last junction
      self1.start1 = self1.overs[0][0] == 0
      self1.start2 = self1.overs[0][1] == 0
      self1.end1 = self1.overs[-1][0] == len(self1.tx_obj1.exons)-1
      self1.end2 = self1.overs[-1][1] == len(self1.tx_obj2.exons)-1
      return

    def calculate_overlap(self1):
      """Create the array that describes how junctions overlap"""
      overs = []
      if not self1.tx_obj1.get_range().overlaps(self1.tx_obj2.get_range()): return # if they dont overlap wont find anything
      for i in range(0,len(self1.tx_obj1.exons)):
        for j in range(0,len(self1.tx_obj2.exons)):
          osize = self1.tx_obj1.exons[i].rng.overlap_size(self1.tx_obj2.exons[j].rng)
          ofrac = 0
          if osize > 0:
            ofrac = min(float(osize)/float(self1.tx_obj1.exons[i].rng.length())\
                       ,float(osize)/float(self1.tx_obj2.exons[j].rng.length()))

          if self1.tx_obj1.get_exon_count() == 1 or self1.tx_obj2.get_exon_count == 1:
            # use single exon rules
            if osize >= self1.single_minover and ofrac >= self1.single_frac:
              #print 'single exon match'
              overs.append([i,j])
          else: # for multi exons
            if i == 0 or j == 0 or i == len(self1.tx_obj1.exons)-1 or j == len(self1.tx_obj2.exons)-1:
              #its on an end
              if osize >= self1.multi_minover and ofrac >= self1.multi_endfrac:
                #print 'end exon match'
                overs.append([i,j])
            #else its a middle
            elif osize >= self1.multi_minover and ofrac >= self1.multi_midfrac:
              #print 'mid exon match'
              overs.append([i,j])
      #print overs
      self1.overs = overs

class JunctionOverlap:
    """Class for describing the overlap of junctions between transcripts
 
    This should probably be not a child.

    :param tx_obj1: transcript1
    :param tx_obj2: transcript2
    :param tolerance: how far before its no longer a matched junction
    :type tx_obj1: Transcript
    :type tx_obj2: Transcript
    :type tolerance: int
    """
    def __init__(self1,tx_obj1,tx_obj2,tolerance=0):
      self1.tx_obj1 = tx_obj1
      self1.tx_obj2 = tx_obj2
      self1.tolerance = tolerance
      self1.overs = [] # gets set by calculate_overlap()
      self1.calculate_overlap()
      if len(self1.overs) == 0: return None# nothing to analyze
      self1.analyze_overs()

    def __nonzero__(self1):
      if len(self1.overs) > 0: return True
      return False

    def match_junction_count(self1):
      return len(self1.overs)

    def is_subset(self1):
      """Return value if tx_obj2 is a complete subset of tx_obj1 or tx_obj1 is a complete subset of tx_obj2

      values:

      * Return 1: Full overlap (mutual subests) 
      * Return 2: two is a subset of one
      * Return 3: one is a subset of two
      * Return False if neither is a subset of the other
      """
      if len(self1.overs) == 0: return False
      if len(self1.dif1) > 0: # make sure they are consecutive if more than one
        if max(self1.dif1) != 1 or max(self1.dif2) != 1: return False
      onecov = self1.start1 and self1.end1
      twocov = self1.start2 and self1.end2
      if onecov and twocov:
        return 1
      elif twocov: return 2
      elif onecov: return 3
      return False

    def is_full_overlap(self1):
      """True if its a full overlap

      :return: True if its a full overlap
      :rtype: bool
      """
      if len(self1.overs) == 0: return False
      if len(self1.dif1) > 0:
        if max(self1.dif1) != 1 or max(self1.dif2) != 1: return False
      if self1.start1 and self1.end1 and self1.start2 and self1.end2:
        return True
      return False

    def is_compatible(self1):
      """Return True if the transcripts can be combined together

      :return: True if we can combine
      :rtype: bool
      """
      if len(self1.overs) == 0: return False
      if len(self1.dif1) > 0:
        if max(self1.dif1) != 1 or max(self1.dif2) != 1: return False
      # If we are still here it is a single run
      if (self1.start1 or self1.start2) and (self1.end1 or self1.end2):
        return True
      return False

    def analyze_overs(self1):
      """A helper function to prepare values describing overlaps"""
      #check for full overlap first
      self1.dif1 = [self1.overs[i][0]-self1.overs[i-1][0] for i in range(1,len(self1.overs))]
      self1.dif2 = [self1.overs[i][1]-self1.overs[i-1][1] for i in range(1,len(self1.overs))]
      #see if it starts and ends on first or last junction
      self1.start1 = self1.overs[0][0] == 0
      self1.start2 = self1.overs[0][1] == 0
      self1.end1 = self1.overs[-1][0] == len(self1.tx_obj1.junctions)-1
      self1.end2 = self1.overs[-1][1] == len(self1.tx_obj2.junctions)-1
      return

    def calculate_overlap(self1):
      """Create the array that describes how junctions overlap"""
      overs = []
      if not self1.tx_obj1.get_range().overlaps(self1.tx_obj2.get_range()): return # if they dont overlap wont find anything
      for i in range(0,len(self1.tx_obj1.junctions)):
        for j in range(0,len(self1.tx_obj2.junctions)):
          if self1.tx_obj1.junctions[i].overlaps(self1.tx_obj2.junctions[j],self1.tolerance):
            overs.append([i,j])
      self1.overs = overs

class Junction:
  """ class to describe a junction

  :param rng_left: left side of junction
  :param rng_right: right side of junction
  :type rng_left: GenomicRange
  :type rng_right: GenomicRange
  """
  def __init__(self,rng_left=None,rng_right=None):
    self.left = rng_left
    self.right = rng_right
    self.left_exon = None
    self.right_exon = None
  def dump_serialized(self):
    """Get string representation of the junction

    :return: serialized object
    :rtype: string
    """
    return pickle.dumps(self)
  def load_serialized(self,instr):
    """load the string

    :param instr:
    :type instr: string
    """
    self = pickle.loads(instr)
  def get_string(self):
    """A string representation of the junction

    :return: string represnetation
    :rtype: string
    """
    return self.left.chr+':'+str(self.left.end)+'-'+self.right.chr+':'+str(self.right.start)
  def get_left_exon(self):
    """ get the exon to the left of the junction

    :return: left exon
    :rtype: Exon
    """
    return self.left_exon
  def get_right_exon(self):
    """ get the exon to the right of the junction

    :return: right exon
    :rtype: Exon or GenomicRange
    """
    return self.right_exon
  def get_range_string(self):
    """Another string representation of the junction.  these may be redundant."""
    return self.left.chr+":"+str(self.left.end)+'/'+self.right.chr+":"+str(self.right.start)
  def set_left(self,rng):
    """ Assign the leftmost range"""
    self.left = rng
  def set_right(self,rng):
    """ Assign the right most range"""
    self.right = rng
  def equals(self,junc):
    """test equality with another junction"""
    if self.left.equals(junc.left): return False
    if self.right.equals(junc.right): return False
    return True
  def overlaps(self,junc,tolerance=0):
    """see if junction overlaps with tolerance""" 
    if not self.left.overlaps_with_padding(junc.left,tolerance): return False
    if not self.right.overlaps_with_padding(junc.right,tolerance): return False
    return True
  def cmp(self,junc,tolerance=0):
    """ output comparison and allow for tolerance if desired

    * -1 if junc comes before self
    * 1 if junc comes after self
    * 0 if overlaps
    * 2 if else

    :param junc:
    :param tolerance: optional search space (default=0, no tolerance)
    :type junc: Junction
    :type tolerance: int
    :return: value of comparison
    :rtype: int
    """
    if self.overlaps(junc,tolerance):
      return 0 #equal
    if self.left.chr == junc.right.chr:
      if self.left.start > junc.right.start:
        return -1 #comes before
    if self.right.chr == junc.left.chr:
      if self.right.start < junc.right.start:
        return 1 #comes after
    return 2

  def set_exon_left(self,ex):
    """assign the left exon"""
    self.left_exon = ex
    ex.right_junc = self
  def set_exon_right(self,ex):
    """assign the right exon"""
    self.right_exon = ex
    ex.left_junc = self

class Exon:
  """class to describe an exon

  :param rng:
  :type rng: GenomicRange
  """
  def __init__(self,rng=None):
    self.rng = rng
    self.left_junc = None
    self.right_junc = None
    self._is_leftmost = False #bool is it a start or end
    self._is_rightmost = False
  def dump_serialized(self):
    return pickle.dumps(self)
  def load_serialized(self,instr):
    self = pickle.loads(instr)
  def get_range(self):
    return self.rng
  def get_length(self):
    return self.rng.length()
  def set_left_junc(self,junc):
    self.left_junc = junc
    junc.set_right_exon = self
  def set_right_junc(self,junc):
    self.right_junc = junc
    junc.set_left_exon = self
  def set_is_leftmost(self,boo=True): 
    self._is_leftmost = boo # is leftmost
  def set_is_rightmost(self,boo=True): 
    self._is_rightmost = boo   #is rightmost

def _mode(mylist):
  counts = [mylist.count(x) for x in mylist]
  maxcount = max(counts)
  avg = sum([float(x) for x in mylist])/len(mylist)
  #print counts 
  dist = [abs(float(x)-avg) for x in mylist]
  best_list = []
  best_dist = []
  for i in range(0,len(mylist)):
    counts[i] == maxcount
    best_list.append(mylist[i])
    best_dist.append(dist[i])
  abs_best_dist = min(best_dist)
  for i in range(0,len(best_dist)):
    if best_dist[i] == abs_best_dist: 
      return best_list[i]
  sys.stderr.write("Warning: trouble finding best\n")
  return best_list[0]


def trim_ordered_range_list(ranges,start,finish):
  """A function to help with slicing a mapping
     Start with a list of ranges and get another list of ranges constrained by start (0-indexed) and finish (1-indexed)

     :param ranges: ordered non-overlapping ranges on the same chromosome
     :param start: start 0-indexed
     :param finish: ending 1-indexed
     :type ranges: GenomicRange []
     :type start: Int
     :type finish: Int
     :return: non-overlapping ranges on same chromosome constrained by start and finish
     :rtype: GenomicRange []
  """
  z = 0
  keep_ranges = []
  for inrng in self.ranges:
    z+=1
    original_rng = inrng
    rng = inrng.copy() # we will be passing it along and possibly be cutting it
    done = False;
    #print 'exon length '+str(rng.length())
    if start >= index and start < index+original_rng.length(): # we are in this one
      rng.start = original_rng.start+(start-index) # fix the start
      #print 'fixstart '+str(original_rng.start)+' to '+str(rng.start)
    if finish > index and finish <= index+original_rng.length():
      rng.end = original_rng.start+(finish-index)-1
      done = True
      #print 'fixend '+str(original_rng.end)+' to '+str(rng.end)
 
    if finish <= index+original_rng.length(): # we are in the last exon we need
      index+= original_rng.length()
      keep_ranges.append(rng)
      break
    if index+original_rng.length() < start: # we don't need any bases from this
      index += original_rng.length()
      continue # we don't use this exon
    keep_ranges.append(rng)
    index += original_rng.length()
    if index > finish: break
    if done: break
  return keep_ranges

