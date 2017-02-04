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

  def _shift_index(self,val):
    """return a value to what its start should be given a index-1 value"""
    return val
  
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
      nrng = type(self)(self.chr,self._shift_index(range2.end+1),self.end,self._options)
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

  def _shift_index(self,val):
    """return a value to what its start should be given a index-1 value"""
    return val-1

  def copy(self):
    """Override copy to make another copy

    :returns: a new copy of this object
    :rtype: Bed

    """
    return type(self)(self.chr,self.start-1,self.end,self._options)

  def __str__(self):
    return str(self.chr)+"\t"+str(self.start)+"\t"+str(self.end)+"\t"+str(self._options)

class Locus:
  """A Locus is a colloction of GenomicRanges that fall 
     within some distance of one another"""
  def __init__(self):
    self.range = None
    self.members = []
    self.use_direction = False #If we want to use direction, input ranges must have the direction set
  def set_use_direction(self,inbool):
    """Set to true if you want all locus members to share the same direction

    :param inbool:
    :type inbool: bool

    """
    self.use_direction = inbool
  def add_member(self,grange):
    """Add a genomic range to the locus

    :param grange:
    :type grange: GenomicRange

    """
    if self.use_direction and not grange.direction:
      sys.stderr.write("ERROR if using direction then direction of input members must be set\n")
      sys.exit()
    # Get range set properly
    if not self.range:
      self.range = GenomicRange(grange.chr,grange.start,grange.end)
      if self.use_direction:
        self.range.set_direction(grange.get_direction())
    elif self.range.chr != grange.chr:
      sys.stderr.write("WARNING cannot add member with chromosomes are not equal\n")
      return False
    elif self.use_direction and self.range.direction != grange.direction:
      sys.stderr.write("WARNING cannot add member with different directions\n")
      return False
    else:
      if grange.start < self.range.start:  self.range.start = grange.start
      if grange.end > self.range.end: self.range.end = grange.end
    self.members.append(grange)

class Loci:
  """multiple locus combined together when new members are added
     based on parameters"""
  def __init__(self):
    self.loci = []
    self.overhang = 0
    self.use_direction = False
    self.verbose = False
    return
  def set_minimum_distance(self,over):
    """In preparation for combining loci specify how many basepairs they may be separated by but still get  merged"""
    self.overhang = over
  def set_use_direction(self,inbool):
    """ Do we want to only combine loci when they have the same direction, if so, set to True"""
    self.use_direction = inbool
  def add_locus(self,inlocus):
    """ Adds a locus to our loci, but does not go through an update our locus sets yet"""
    if self.use_direction == True and inlocus.use_direction == False:
      sys.stderr.write("ERROR if using the direction in Loci, then every locus added needs use_direction to be True\n")
      sys.exit()
    self.loci.append(inlocus)
    return
  def update_loci(self):
    """Goes through and combines loci until we have one set meeting our overlap definition"""
    # Create sub-loci for each chromosome
    lbc = {}
    chroms = sorted([x.range.chr for x in self.loci])
    for chrom in chroms: lbc[chrom] = Loci()
    for x in self.loci: lbc[x.range.chr].add_locus(x)
    for chrom in sorted(lbc.keys()):
      if self.verbose: 
        lbc[chrom].verbose = True
        sys.stderr.write(chrom+"\n")
      lbc[chrom].overhang = self.overhang
      lbc[chrom].use_direction = self.use_direction
      lbc[chrom].merge_down_loci()
    self.loci = []
    for chrom in sorted(lbc.keys()):
      for locus in lbc[chrom].loci:  self.loci.append(locus)
  def merge_down_loci(self):
    """Called internally to make loci overlapping into one set"""
    old_locus_size = -1
    z = 0
    while len(self.loci) != old_locus_size:
      z+=1
      old_locus_size = len(self.loci)
      locus_size = len(self.loci)
      if self.verbose:
        sys.stderr.write(str(locus_size)+" Combining down loci step "+str(z)+"       \r")
      combined = set()
      for i in range(0,locus_size):
        if i in combined: continue
        for j in range(i+1,locus_size):
          if self.loci[i].range.overlaps_with_padding(self.loci[j].range,self.overhang):
            if self.use_direction and self.loci[i].range.direction != self.loci[j].range.direction:  continue
            for obj in self.loci[j].members:
              self.loci[i].add_member(obj)
            combined.add(j)
            break
      newloci = []
      for i in range(0,locus_size):
        if i not in combined:
          newloci.append(self.loci[i])
      self.loci = newloci
    if self.verbose:
      sys.stderr.write("Finished combining down "+str(len(self.loci))+" loci in "+str(z)+" steps   \n")
    return

def sort_ranges(inranges):
  """from an array of ranges, make a sorted array of ranges

  :param inranges: List of GenomicRange data
  :type inranges: GenomicRange[]
  :returns: a new sorted GenomicRange list
  :rtype: GenomicRange[]

  """
  return sorted(inranges,key=lambda x: (x.chr,x.start,x.end,x.direction))
  
def merge_ranges(inranges,already_sorted=False):
  """from a list of genomic range or bed entries, whether or not they are already sorted, 
     make a flattend range list of ranges where if they overlapped, they are now joined
     (not yet) The new range payloads will be the previous ranges

  :param inranges:
  :param already_sorted: has this already been sorted (defaults to False)
  :type inranges: GenomicRange[]
  :type already_sorted: bool

  :return: sorted ranges
  :rtype: GenomicRange[]

  """
  if not already_sorted: inranges = sort_ranges(inranges)
  prev = None
  outputs = []
  merged = False
  for rng in inranges:
    #nrng = rng.copy()
    #nrng.set_payload([])
    #nrng.get_payload().append(rng)
    merged = False
    if len(outputs) > 0:
      if rng.overlaps(outputs[-1]) or rng.adjacent(outputs[-1]):
        nrng = rng.merge(outputs[-1])
        #nrng.set_payload(prev.get_payload())
        #nrng.get_payload().append(rng)
        outputs[-1] = nrng
        merged = True
    if not merged:
      outputs.append(rng.copy())
    #prev = nrng
  #if not merged: outputs.append(prev)
  return sort_ranges(outputs)

def pad_ranges(inranges,padding,chr_ranges=None):
  """Add the specfied amount onto the edges the transcripts

  :param inranges: List of genomic ranges in Bed o GenomicRange format.
  :param padding: how much to add on
  :param chr_ranges: looks like the list of ranges within which to pad
  :type inranges: GenomicRange[]
  :type padding: int
  :type chr_ranges: 

  """
  if not inranges: return
  outranges = []
  if len(inranges) == 0: return outranges
  chr = {}
  if chr_ranges:
    for b in chr_ranges:
      chr[b.chr] = b
  for rng in inranges:
    newstart = rng.start - padding
    newend = rng.end + padding
    if rng.chr in chr:
      if newstart < chr[rng.chr].start: newstart = chr[rng.chr].start
      if newend > chr[rng.chr].end: endstart = chr[rng.chr].end
    nrng = rng.copy()
    nrng.start = newstart
    nrng.end = newend
    outranges.append(nrng)
  return sort_ranges(outranges)

def subtract_ranges(r1s,r2s,already_sorted=False):
  """Subtract multiple ranges from a list of ranges

  :param r1s: range list 1
  :param r2s: range list 2
  :param already_sorted: default (False)
  :type r1s: GenomicRange[]
  :type r2s: GenomicRange[]

  :return: new range r1s minus r2s
  :rtype: GenomicRange[]
  """
  from seqtools.stream import MultiLocusStream
  if not already_sorted:
    r1s = merge_ranges(r1s)
    r2s = merge_ranges(r2s)
  outputs = []
  mls = MultiLocusStream([BedArrayStream(r1s),BedArrayStream(r2s)])
  tot1 = 0
  tot2 = 0
  for loc in mls:
    #[beds1,beds2] = loc.get_payload()
    v = loc.get_payload()
    #print v
    [beds1,beds2] =v
    beds1 = beds1[:]
    beds2 = beds2[:]
    if len(beds1)==0:
      continue
    if len(beds2)==0:
      outputs += beds1
      continue
    #this loop could be made much more efficient
    mapping = {} #keyed by beds1 index stores list of overlaping beds2 indecies
    for i in range(0,len(beds1)):
      mapping[i] = []
    beds2min = 0
    beds2max = len(beds2)
    for i in range(0,len(beds1)):
      for j in range(beds2min,beds2max):
        cmpval = beds1[i].cmp(beds2[j])
        if cmpval == -1:
          beds2min = j+1
        elif cmpval == 0:
          mapping[i].append(j)
        else:
          break
    for i in range(0,len(beds1)):
      if len(mapping[i])==0: outputs += beds1
      else:
        outputs += subtract_range_array(beds1[i],[beds2[j] for j in mapping[i]],is_sorted=True)
    #while len(beds2) > 0:
    #  b2 = beds2.pop(0)
    #  vs = [x.subtract(b2) for x in beds1]
    #  tot = []
    #  for res in vs:
    #    tot = tot + res
    #  beds1 = tot
    #print "subtract "+str(len(beds1))+"\t"+str(len(beds2))
    #print beds1[0].get_range_string()
  #outputs = merge_ranges(outputs)
  #print [x.get_range_string() for x in outputs]

  return merge_ranges(outputs)

def union_range_array(bed1,beds2,payload=None,is_sorted=False):
  """ Does not do a merge if the payload has been set

  :param bed1:
  :param bed2:
  :param payload: payload=1 return the payload of bed1 on each of the union set, payload=2 return the payload of bed2 on each of the union set, payload=3 return the payload of bed1 and bed2 on each of the union set
  :param is_sorted:
  :type bed1: GenomicRange
  :type bed2: GenomicRange
  :type payload: int
  :type is_sorted: bool
  """
  if not is_sorted: beds2 = sort_ranges(beds2)
  output = []
  for bed2 in beds2:
    cval = bed2.cmp(bed1)
    #print str(cval)+" "+bed1.get_range_string()+" "+bed2.get_range_string()
    if cval == -1: continue
    elif cval == 0:
      output.append(bed1.union(bed2))
      if payload==1:
        output[-1].set_payload(bed1.get_payload())
      if payload==2:
        output[-1].set_payload(bed2.get_payload())
    elif cval == 1: break
  if payload: return sort_ranges(output)
  return merge_ranges(output)

def subtract_range_array(bed1,beds2,is_sorted=False):
  """subtract several ranges from a range, returns array1 - (all of array2)

  :param bed1: 
  :param beds2: subtract all these beds from bed1
  :param is_sorted: has it been sorted already? Default (False)
  :type bed1: Bed or GenomicRange
  :type beds2: Bed[] or GenomicRange[]
  :param is_sorted: bool

  """
  if not is_sorted: beds2 = sort_ranges(beds2)
  output = [bed1.copy()]  
  mink = 0
  for j in range(0,len(beds2)):
    temp = []
    if mink > 0: temp = output[0:mink]
    for k in range(mink,len(output)):
      cmpv = output[k].cmp(beds2[j])
      if cmpv ==-1: mink=k
      temp += output[k].subtract(beds2[j])
    #for nval in [x.subtract(beds2[j]) for x in output]:
    #  temp += nval
    output = temp
  return output

def string_to_genomic_range(rstring):
  """ Convert a string to a genomic range

  :param rstring: string representing a genomic range chr1:801-900
  :type rstring:
  :returns: object representing the string
  :rtype: GenomicRange
  """
  m = re.match('([^:]+):(\d+)-(\d+)',rstring)
  if not m: 
    sys.stderr.write("ERROR: problem with range string "+rstring+"\n")
  return GenomicRange(m.group(1),int(m.group(2)),int(m.group(3)))

def sort_genomic_ranges(rngs):
  """sort multiple ranges"""
  return sorted(rngs, key=lambda x: (x.chr, x.start, x.end))

def ranges_to_coverage(rngs,threads=1):
  """take a list of ranges as an input
  output a list of ranges and the coverage at each range
  :param rngs: bed ranges on a single chromosome. not certain about that single chromosome requirement
  :type rngs: GenomicRange[] or Bed[]
  :param threads: Not currently being used
  :type threads: int

  :return: out is the non-overlapping bed ranges with the edition of depth
  :rtype: GenomicRange[]
  """
  def do_chr(rngs):
    """do one chromosomes sorting
    :param rngs:
    :type rngs: GenomicRange[]
    """
    #starts = sorted(range(0,len(rngs)), key=lambda x: rngs[x].start)
    #print starts
    #ends = sorted(range(0,len(rngs)), key=lambda x: rngs[x].end)
    start_events = [x.start for x in rngs]
    end_events = [x.end+1 for x in rngs]
    indexed_events = {}
    for e in start_events:
      if e not in indexed_events: indexed_events[e] = {'starts':0,'ends':0}
      indexed_events[e]['starts']+=1
    for e in end_events:
      if e not in indexed_events: indexed_events[e] = {'starts':0,'ends':0}
      indexed_events[e]['ends']+=1
    cdepth = 0
    pstart = None
    pend = None
    outputs = []
    ordered_events = sorted(indexed_events.keys())
    for loc in ordered_events:
      prev_depth = cdepth # where we were
      # see where we are before the change
      cdepth += indexed_events[loc]['starts']
      cdepth -= indexed_events[loc]['ends']
      if prev_depth > 0 and prev_depth != cdepth:
        outputs.append([rngs[0].chr,pstart,loc-1,prev_depth]) # output what was before this if we are in something
      if prev_depth != cdepth or cdepth == 0:
        pstart = loc
    #print outputs
    return outputs
  
  class Queue:
    """Simple class to be able to use get function to retreive a value"""
    def __init__(self,val):
      self.val = [val]
    def get(self):
      return self.val.pop(0)

  ### START MAIN ####
  srngs = sort_genomic_ranges(rngs)
  # get the leftmost unique range
  chr = srngs[0].chr
  buffer = []
  results = []
  prelim = []
  for b in srngs:
    if b.chr != chr:
      rs = do_chr(buffer[:])
      for r in rs:  
        results.append(GenomicRange(r[0],r[1],r[2]))
        results[-1].set_payload(r[3])
      buffer = []
    buffer.append(b)
    chr = b.chr
  if len(buffer) > 0:
    rs = do_chr(buffer[:])
    for r in rs: 
      results.append(GenomicRange(r[0],r[1],r[2]))
      results[-1].set_payload(r[3])
  return results

class BedArrayStream:
  """Make a stream from a bedarray

  Read as an interator or with read_entry()

  :param bedarray:
  :type bedarray: Bed[]
  """
  def __init__(self,bedarray):
    self.prev = None
    self.curr_ind = 0
    self.bedarray = bedarray
  def read_entry(self):
    """get the next value from the array, and set internal iterator so next call will be next entry

    :return: The next GenomicRange entry
    :rtype: GenomicRange
    """
    if len(self.bedarray) <= self.curr_ind: return None
    val = self.bedarray[self.curr_ind]
    self.curr_ind += 1
    return val
  def next(self):
    """ call read_entry() from inside this iterator"""
    r = self.read_entry()
    if not r: raise StopIteration
    else:
      return r
  def __iter__(self):
    return self

class BedStream:
  """Make a stream from a handle, keep it as an iterator

  :param fh: readable file handle or stream
  :type fh: handle
  """
  def __init__(self,fh):
    self.fh = fh
  def read_entry(self):
    """read the next bed entry from the stream"""
    line = self.fh.readline()
    if not line: return None
    m = re.match('([^\t]+)\t(\d+)\t(\d+)\t*(.*)',line.rstrip())
    if not m:
      sys.stderr.write("ERROR: unknown line in bed format file\n"+line+"\n")
      sys.exit()
    g = GenomicRange(m.group(1),int(m.group(2))+1,int(m.group(3)))
    if len(m.group(4)) > 0:
      g.set_payload(m.group(4))
    return g
  def next(self):
    r = self.read_entry()
    if not r: raise StopIteration
    else:
      return r
  def __iter__(self):
    return self

