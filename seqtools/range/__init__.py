""" These classes are to help deal with ranged data
    things associated with those coordinates. """
import sys, re
from collections import namedtuple


class RangeGeneric(object):
  """A generic range object

     Use 1-index start and 1-index end

     Slicing is permitted and returns a new RangeGeneric but keep in mind
     the slicing is 1-index as well for this object
  """
  def __init__(self,start,end,options=None):
    if not options: options = RangeGeneric.Options()
    self._start = start
    self._end = end
    self._options = options

  @staticmethod
  def Options(**kwargs):
     """Create a new options namedtuple with only allowed keyword arguments"""
     attributes = ['payload']
     Opts = namedtuple('Opts',attributes)
     if not kwargs: return Opts(**dict([(x,None) for x in attributes]))
     return Opts(**dict(kwargs))
  
  def __len__(self):
    return self._end-self._start+1

  @property
  def start(self):
    return self._start

  @property
  def start0(self):
    return self._start-1

  @property
  def start1(self):
    return self._start

  @property
  def end(self):
    return self._end

  def copy(self):
    return type(self)(self._start,self._end,self._options)

  @property
  def payload(self):
    """Returns the payload, whatever it may be"""
    return self._options.payload

  def set_payload(self,inpay):
    """Set the payload.  Stored in a list to try to keep it as a reference

    :param inpay: payload input
    """
    self._options = options._replace(payload=inpay)

  def __iter__(self):
    """Lets try to do it by makinga n iterator"""
    for i in range(self._start,self._end+1):
      yield i

  def __getitem__(self,key):
    if key.step:
      raise ValueError('a step should not be used when slicing a sequence range')
    if isinstance(key,slice):
      if key.start > self._end: return None
      if key.stop < self._start: return None
      return RangeGeneric(max(key.start,self._start),min(key.stop,self._end),self._options)

  def __setitem__(self,key):
    return
  def __delitem__(self,key):
    return
  def __str__(self):
    return str(self._start)+'-'+str(self._end)+' '+str(self._options)


class GenomicRange(RangeGeneric):
  """A basic class for keeping genomic range data.  It is 1-indexed for both start and end.

  It can carry directional information, but this information is 
  not required for equality and adjacency etc.

  :param chr: chromosome name
  :param start: 1-indexed starting base
  :param end: 1-indexed ending base
  :param options: namedtuple of parameters
  :param options.dir: optional direction
  :param options.payload: optional payload
  :type chr: char
  :type start: int
  :type end: int
  :type options: namedtuple()
  :type options.dir: Char
  :type options.payload: Object

  """
  def __init__(self,chr,start,end,options=None):
    if not options: 
      options = GenomicRange.Options()
    super(GenomicRange,self).__init__(start,end,options)
    self._chr = chr

  def __getitem__(self,key):
    if key.step:
      raise ValueError('a step should not be used when slicing a sequence range')
    if isinstance(key,slice):
      if key.start > self._end: return None
      if key.stop < self._start: return None
      return GenomicRange(self._chr,max(key.start,self.start1),min(key.stop,self._end),self._options)

  @staticmethod
  def Options(**kwargs):
     """Create a new options namedtuple with only allowed keyword arguments"""
     attributes = ['payload','dir']
     Opts = namedtuple('Opts',attributes)
     if not kwargs: return Opts(**dict([(x,None) for x in attributes]))
     return Opts(**dict(kwargs))

  @property
  def chr(self):
    return self._chr

  def __str__(self):
    return self.get_range_string()+' '+str(self._options)

  def copy(self):
    """Create a new copy of selfe.  does not do a deep copy for payload

    :return: copied range
    :rtype: GenomicRange

    """
    return type(self)(self.chr,self.start,self.end,self._options)
    #return GenomicRange(self._chr,self._start,self._end,self._options)

  def get_range(self):
    """For compatability with some range-based tools that need to call this function

    :return: this object
    :rtype: GenomicRange
    """
    return self

  def get_bed_array(self):
    """Return a basic three meber bed array representation of this range

    :return: list of [chr,start (0-indexed), end (1-indexed]
    :rtype: list
    """
    arr = [self._chr,self._start-1,self._end]
    if self._options.dir:
      arr.append(self._options.dir)
    return arr 

  @property
  def direction(self):
    """return the direction

    :return: the direction or strand +/- (or None if not set)
    :rtype: char

    """
    return self._options.dir

  def set_direction(self,dir):
    """ set he direction

    :param dir: direction + or -
    :type dir: char
    """
    self._options = self._options._replace(dir=dir)

  def equals(self,gr):
    """ check for equality. does not consider direction

    :param gr: another genomic range
    :type gr: GenomicRange
    :return: true if they are the same, false if they are not
    :rtype: bool
    """
    if self.chr == gr.chr and self.start1 == gr.start1 and self.end == gr.end:
      return True
    return False

  def get_range_string(self):
    """ get the range string represetation. similar to the default input for UCSC genome browser

    :return: representation by string like chr2:801-900
    :rtype: string
    """
    return self.chr+":"+str(self.start1)+"-"+str(self.end)

  def get_bed_coordinates(self):
    """ Same as get bed array.
        These are the 0-indexed start, 1-indexted stop coordinates

    :return: bed array [chr,start-1, end]
    """
    return [self.chr,self.start0,self.end]

  def get_genomic_coordinates(self):
    """These are the 1-indexed coordiantes in list form

    :return: list of coords [chr, start (1-indexed), end(1-indexed)
    :rtype: list

    """
    return [self.chr,self.start1,self.end]

  def adjacent(self,rng2):
    """ Test for adjacency.  

    :param rng2:
    :param use_direction: false by default
    :param type: GenomicRange
    :param type: use_direction
    """
    if self.chr != rng2.chr: return False
    if self.direction != rng2.direction and use_direction: return False
    if self.end == rng2.start1-1:  return True
    if self.start1-1 == rng2.end: return True
    return False

  def overlaps(self,in_genomic_range,padding=0):
    """do the ranges overlap?

    :param in_genomic_range: range to compare to
    :param padding: add to the ends this many (default 0)
    :type in_genomic_range: GenomicRange
    :type padding: int

    :return: True if they overlap
    :rtype: bool

    """
    if padding > 0:
      in_genomic_range = GenomicRange(in_genomic_range.chr,max([1,in_genomic_range.start1-padding]),in_genomic_range.end+padding)
    if self.chr != in_genomic_range.chr:
      return False
    if self.end < in_genomic_range.start1:
      return False
    if in_genomic_range.end < self.start1:
      return False
    if self.start1 > in_genomic_range.end:
      return False
    if in_genomic_range.start1 > self.end:
      return False
    if self.start1 <= in_genomic_range.start1 and self.end >= in_genomic_range.start1:
      return True
    if self.start1 <= in_genomic_range.end and self.end >= in_genomic_range.end:
      return True
    if self.start1 >= in_genomic_range.start1 and self.end <= in_genomic_range.end:
      return True
    if self.start1 <= in_genomic_range.start1 and self.end >= in_genomic_range.end:
      return True
    if in_genomic_range.start1 <= self.start1 and in_genomic_range.end >= self.start1:
      return True
    if in_genomic_range.start1 <= self.end and in_genomic_range.end >= self.end:
      return True
    sys.stderr.write("overlaps: unprogrammed error\n")
    return False

  def overlap_size(self,in_genomic_range):
    """ The size of the overlap

    :param in_genomic_range: the range to intersect
    :type in_genomic_range: GenomicRange

    :return: count of overlapping bases
    :rtype: int
    """
    if self.chr != in_genomic_range.chr:
      return 0
    if self.end < in_genomic_range.start1:
      return 0
    if in_genomic_range.end < self.start1:
      return 0
    if self.start1 > in_genomic_range.end:
      return 0
    if self.start1 >= in_genomic_range.start1 and self.end <= in_genomic_range.end:
      return self.end-self.start1+1
    if self.start1 <= in_genomic_range.start1 and self.end >= in_genomic_range.end:
      return in_genomic_range.end-in_genomic_range.start1+1
    if self.start1 <= in_genomic_range.start1 and self.end >= in_genomic_range.start1:
      return self.end-in_genomic_range.start1+1
    if self.start1 <= in_genomic_range.end and self.end >= in_genomic_range.end:
      return in_genomic_range.end-self.start1+1
    if in_genomic_range.start1 <= self.start1 and in_genomic_range.end >= self.start1:
      return in_genomic_range.end-self.start1+1
    if in_genomic_range.start1 <= self.end and in_genomic_range.end >= self.end:
      return self.end-in_genomic_range.start1+1
    sys.stderr.write("overlap_size: unprogrammed error\n")
    return 0

  def merge(self,range2): 
    """merge this bed with another bed to make a longer bed.  Returns None if on different chromosomes.

    keeps the options of this class (not range2)

    :param range2:
    :type range2: GenomicRange

    :return: bigger range with both
    :rtype: GenomicRange

    """
    if self.chr != range2.chr:
      return None
    o = type(self)(self.chr,min(self.start,range2.start),max(self.end,range2.end),self._options)
    return o

  def intersect(self,range2):
    """Return the chunk they overlap as a range.

    options is passed to result from this object

    :param range2:
    :type range2: GenomicRange

    :return: Range with the intersecting segement, or None if not overlapping
    :rtype: GenomicRange

    """
    if not self.overlaps(range2): return None
    return type(self)(self.chr,max(self.start,range2.start),min(self.end,range2.end),self._options)

  def cmp(self,range2,overlap_size=0):
    """the comparitor for ranges

     * return 1 if greater than range2
     * return -1 if less than range2
     * return 0 if overlapped

    :param range2:
    :param overlap_size: allow some padding for an 'equal' comparison (default 0)
    :type range2: GenomicRange
    :type overlap_size: int

    """
    if self.overlaps(range2,padding=overlap_size): return 0
    if self.chr < range2.chr: return -1
    elif self.chr > range2.chr: return 1
    elif self.end < range2.start1: return -1
    elif self.start1 > range2.end: return 1
    sys.stderr.write("ERROR: cmp function unexpcted state\n")
    sys.exit()
    return 0

  def subtract(self,range2):
    """Take another range, and list of ranges after removing range2, keep options from self

    :param range2:
    :type range2: GenomicRange
    :return: List of Genomic Ranges
    :rtype: GenomicRange[]

    """
    outranges = []
    if self.chr != range2.chr:
      outranges.append(self.copy())
      return outranges
    if not self.overlaps(range2):
      outranges.append(self.copy())
      return outranges
    if range2.start1 <= self.start1 and range2.end >= self.end:
      return outranges #delete all
    if range2.start1 > self.start1: #left side
      nrng = type(self)(self.chr,self.start,range2.start1-1,self._options)
      outranges.append(nrng)
    if range2.end < self.end: #right side
      #ugly addon to make it work for either 0 or 1 index start
      nrng = type(self)(self.chr,range2.end+1-(self.start1-self.start),self.end,self._options)
      outranges.append(nrng)
    return outranges

  def equals(self,rng):
    if self.chr != rng.chr: return False
    if self.start1 != rng.start1: return False
    if self.end != rng.end: return False
    return True

  def distance(self,rng):
    """The distance between two ranges.

    :param rng: another range
    :type rng: GenomicRange
    :returns: bases separting, 0 if overlapped or adjacent, -1 if on different chromsomes
    :rtype: int
    """
    if self.chr != rng.chr: return -1
    c = self.cmp(rng)
    if c == 0: return 0
    if c < 0:
      return rng.start1 - self.end-1
    return self.start1 - rng.end-1

def GenomicRangeFromString(range_string,options):
  """Constructor for a GenomicRange object that takes a string"""
  m = re.match('^(.+):(\d+)-(\d+)$',range_string)
  if not m:  
    sys.stderr.write("ERROR bad genomic range string\n"+range_string+"\n")
    sys.exit()
  chr = m.group(1)
  start = int(m.group(2))
  end = int(m.group(3))
  return GenomicRange(chr,start,end,options)

class Bed(GenomicRange):
  """ Bed format is a chromosome, start (0-index), end (1-index). 
      It is a child of GenomicRange but modifies the class 
      to use the 0-based start 1-based end style of a bed file

  :param chrom:
  :param start: 0-indexed
  :param finish: 1-indexed
  :param options: 
  :type chrom: string
  :type start: int
  :type finish: int
  :type options: namedtuple
  """

  def __init__(self,chrom,start,finish,options=None):
    if not options: options = Bed.Options()
    super(Bed,self).__init__(chrom,start+1,finish,options)

  @staticmethod
  def Options(**kwargs):
     """Create a new options namedtuple with only allowed keyword arguments"""
     attributes = ['payload','dir']
     Opts = namedtuple('Opts',attributes)
     if not kwargs: return Opts(**dict([(x,None) for x in attributes]))
     return Opts(**dict(kwargs))

  @property
  def start(self):
    return self._start-1

  def copy(self):
    """Override copy to make another copy

    :returns: a new copy of this object
    :rtype: Bed

    """
    return type(self)(self.chr,self.start,self.end,self._options)

  def __str__(self):
    return str(self.chr)+"\t"+str(self.start)+"\t"+str(self.end)+"\t"+str(self._options)

