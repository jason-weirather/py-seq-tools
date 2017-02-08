""" These classes are to help deal with multiple ranges, i.e.
    Lists of ranges. """
import sys, re
from collections import namedtuple
from seqtools.range import GenomicRange

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
    v = loc.payload
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

def intersect_range_array(bed1,beds2,payload=None,is_sorted=False):
  """ Does not do a merge if the payload has been set

  :param bed1:
  :param bed2:
  :param payload: payload=1 return the payload of bed1 on each of the intersect set, payload=2 return the payload of bed2 on each of the union set, payload=3 return the payload of bed1 and bed2 on each of the union set
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
      output.append(bed1.intersect(bed2))
      if payload==1:
        output[-1].set_payload(bed1.payload)
      if payload==2:
        output[-1].set_payload(bed2.payload)
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

