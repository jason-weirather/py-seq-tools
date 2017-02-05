""" Classes to help stream biological data"""
import sys
from seqtools.range.multi import merge_ranges
from subprocess import Popen, PIPE


class LocusStream:
  """Works for any stream with ordered range bound objects that have the functions.
     LocusStream is a stream itself, and is iterable

   1. read_entry()
   2. get_range()

   Data is not stored as an actual Locus object, but rather in list in the payload of the range covered by the locus

   :param stream: ordered stream with range
   :type stream: Stream
   """

  def __init__(self,stream):
    self.stream = stream
    self.current_range = None
    firstobj = self.stream.read_entry()
    if not firstobj: return
    self.current_range = firstobj.get_range()
    self.current_range.set_payload([firstobj])

  def __iter__(self):
    return self
  def next(self):
    r = self.read_entry()
    if not r: raise StopIteration
    else:
      return r

  def read_entry(self):
    """As long as entires overlap keep putting them together in a list that is
       the payload for a range that describes the bounds of the list

    :return: range with payload list of elements
    :rtype: GenomicRange
    """
    if not self.current_range:
      return None
    output = None
    while True:
      e = self.stream.read_entry()
      if e:
        rng = e.get_range()
        if not rng: continue # continue if nonetype for range
        if rng.overlaps(self.current_range):
          self.current_range.payload.append(e)
          if self.current_range.end < rng.end: self.current_range.end = rng.end
        else: 
          output = self.current_range
          self.current_range = rng
          self.current_range.set_payload([e])
          break
      else:
        output = self.current_range
        self.current_range = None
        break
    return output

class MultiLocusStream:
  """Take an array streams
     Each element should be sorted by position
     Streams need to have this method:

     1. read_entry()
     2. get_range()

     :param streams: list of streams
     :type streams: list
  """
  def __init__(self,streams):
    self.streams = streams
    self.buffers = []
    # seed the buffers
    for i in range(0,len(streams)):
      entry = self.streams[i].read_entry()
      self.buffers.append(entry)
    #self.set_current_range()

  def __iter__(self):
    return self

  def next(self):
    r = self.read_entry()
    if not r: raise StopIteration
    else:
      return r

  def read_entry(self):
    """get the next aggrogate of streams

    :return: range containing a list of entries from each stream that are from the overlapping part
    :rtype: GenomicRange
    """
    output = []
    for i in self.buffers: output.append([])
    rngs = [x.get_range() for x in self.buffers if x]
    #print rngs
    if len(rngs) == 0: return None
    srngs = sorted(rngs,key=lambda x: (x.chr,x.start,x.end))
    mrngs = merge_ranges(srngs)
    current_range = mrngs[0]
    #print current_range.get_range_string()
    got_overlap = True
    while got_overlap == True:
      got_overlap = False
      for i in range(0,len(self.buffers)):
        if not self.buffers[i]: continue #end of this one
        v = self.buffers[i].get_range().cmp(current_range)
        if v==0:
          got_overlap = True
          if self.buffers[i].get_range().overlaps(current_range):
            current_range = current_range.merge(self.buffers[i].get_range())
          output[i].append(self.buffers[i])
          self.buffers[i] = self.streams[i].read_entry()
          #print str([len(x) for x in output])+"\t"+current_range.get_range_string()
        #print str(i)+":"+str(v)
    current_range.set_payload(output)
    return current_range

class GZippedOutputFile:
  """ use gzip utility to compress output

  :param filename: filename to write to
  :type filename: string

  """
  def __init__(self,filename):
    self._fh = open(filename,'w')
    cmd = "gzip"
    self._pipe = Popen(cmd.split(),stdout=self._fh,stdin=PIPE,close_fds=True)
    self._sh = self._pipe.stdin
  def write(self,value):
    self._sh.write(value)
  def close(self):
    self._pipe.communicate()
    self._fh.close()
