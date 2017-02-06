""" Classes to help stream biological data"""
import sys, re, itertools
from seqtools.range.multi import merge_ranges
from subprocess import Popen, PIPE
from format.gpd import GPD
from multiprocessing import Pool, cpu_count
from collections import namedtuple

BufferedLineGeneratorOptions = namedtuple('BufferedStreamOptions',
   ['buffer_size',
    'mapping_function'])

class BufferedLineGenerator:
   """Generic class for a stream that reads a data stream"""
   def __init__(self,stream,options=None):
      if not options: options = BufferedLineGenerator.Options()
      self._options = options
      self._stream = stream
      self._buffer = ''

   @staticmethod
   def Options(**kwargs):
      """ A method for declaring options for the class"""
      construct = BufferedLineGeneratorOptions #IMPORTANT!  Set this
      names = construct._fields
      d = {}
      for name in names: d[name] = None #default values
      """set defaults here"""
      d['buffer_size'] = 1000000
      for k,v in kwargs.iteritems():
         if k in names: d[k] = v
         else: raise ValueError('Error '+k+' is not a property of these options')
      """Create a set of options based on the inputs"""
      return construct(**d)

   def __iter__(self):
      if not self._options.mapping_function:
         return self._gen()
      return itertools.imap(self._options.mapping_function,self._gen())
   def _gen(self):
      while True:
         data = self._buffer
         self._buffer = ''
         read_data = self._stream.read(self._options.buffer_size)
         if read_data: data += read_data
         last = 0
         for m in re.finditer('\n',data):
            yield data[last:m.start()]
            last = m.start()+1
         self._buffer = data[last:]
         if not read_data and last == 0:
            if len(data) > 0:
               yield data
            else:
               break

BufferedLineStreamOptions = BufferedLineGeneratorOptions
"""Just set up an alias for the options constructor since they are same
   between the stream and the generator"""
class BufferedLineStream(object):
   """ This is mostly just a wrapper around the BufferedLineGenerator.
       This class facilitates streaming function calls better.

       Uses same options as buffered line Generator

   """
   def __init__(self,stream,options=None):
      if not options: options = BufferedLineGenerator.Options()
      self._iterator = BufferedLineGenerator(
         stream,
         options=options
      ).__iter__()
      self._previous = None # saving for some future children that can monitor the stream

   @staticmethod
   def Options(**kwargs):
      return BufferedLineGenerator.Options(**kwargs)

   def __iter__(self):
      return self
   def next(self):
      try:
         r = self._iterator.next()
      except StopIteration: r = None
      if not r: raise StopIteration
      else:
         return r

LineObjectStreamOptions = namedtuple('BufferedStreamOptions',
   ['buffer_size'])
class LineObjectStream(BufferedLineStream):
   """This is probably the version of this family of streamers that will
      either see extensive use or get extended

   :param stream: usually a file handle or similar stream
   :param mapping_function: a function that can convert a line into an object of interest
   :param options: a set of parameters to tune the stream like buffer_size
   :type stream: stream or filehandle
   :type mapping_function: function with a one line, single string input, and one object output
   """
   def __init__(self,stream,mapping_function,options=None):
      if not options: options = LineObjectStream.Options()
      super(LineObjectStream,self).__init__(stream,BufferedLineStream.Options(
         buffer_size = options.buffer_size,
         mapping_function = mapping_function
      ))      

   @staticmethod
   def Options(**kwargs):
      """ A method for declaring options for the class"""
      construct = LineObjectStreamOptions #IMPORTANT!  Set this
      names = construct._fields
      d = {}
      for name in names: d[name] = None #default values
      """set defaults here"""
      d['buffer_size'] = 1000000
      for k,v in kwargs.iteritems():
         if k in names: d[k] = v
         else: raise ValueError('Error '+k+' is not a property of these options')
      """Create a set of options based on the inputs"""
      return construct(**d)

class OrderedStream(object):
   """Encapsulates another stream that has the following properties:

      1. It can be traversed with the next() command
      2. The values it returns have a .range property

      It garuntees that you will have an ordered stream

   :param stream:
   :type stream: an iterable stream
   """
   def __init__(self,stream):
      self._stream = stream
      self._previous = None
      self._ln = 0
   def __iter__(self):
      return self
   def next(self):
      r = self._stream.next() # may already be throwing StopIteration
      if not r: raise StopIteration
      rng = r.range
      if not rng: raise ValueError('The objects streamed must have a get_range() function that returns a range')
      if self._previous:
         if self._previous.cmp(rng) > 0:
            raise ValueError('Expected lines to be ordered but they appear not to be ordered on line '+str(self._ln))
      self._previous = rng
      self._ln += 1
      return r

class LocusStream(object):
   """Works for any stream with ordered range bound objects that have the functions.
      LocusStream is a stream itself, and is iterable

    1. next() to access iterator
    2. .range property

    Data is not stored as an actual Locus object, but rather in list in the payload of the range covered by the locus

    :param stream: ordered stream with range
    :type stream: Stream
    """

   def __init__(self,stream):
      self._stream = OrderedStream(stream)
      self._current_range = None
      try:
         firstobj = self._stream.next()
      except StopIteration: firstobj = None
      if not firstobj: return
      self._current_range = firstobj.range
      self._current_range.set_payload([firstobj])

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
      if not self._current_range:
         return None
      output = None
      while True:
         try:
           e = self._stream.next()
         except StopIteration: e = None
         if e:
            rng = e.range
            if not rng: 
               raise ValueError('no range property. it is required in a locus stream') 
            if rng.overlaps(self._current_range):
               self._current_range.payload.append(e)
               if self._current_range.end < rng.end: self._current_range.end = rng.end
            else: 
               output = self._current_range
               self._current_range = rng
               self._current_range.set_payload([e])
               break
         else:
            output = self._current_range
            self._current_range = None
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
