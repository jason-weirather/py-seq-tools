"""Classes to work with sam and bam files"""
import struct, zlib, sys, re, itertools,gzip
from collections import namedtuple
import seqtools.format.sam
from seqtools.format.sam.header import SAMHeader
from seqtools.format.sam.bam import BAM
from subprocess import Popen, PIPE

class BAMFileGeneric(object):
  """Base class for accessing bam files"""
  def __init__(self,filename,options):
    self._options = options
    self._path = filename
    # set up defaults that will be set by children
    self._header_text = ''
    self._n_ref = 0
    self._ref_lengths = {}
    self._ref_names = []
    self._fh = None # must be set by child init after super call

  @property
  def path(self): return self._path

  def _fetch_headers(self):
    """Needs ._fh handle to stream to be set by child"""
    self._header_text, self._n_ref = self._read_top_header()
    self._ref_lengths, self._ref_names = self._read_reference_information()
    self._header = SAMHeader(self._header_text)

  @property
  def header(self):
    return self._header

  def close(self):
    self._fh.close()

  def __iter__(self):
    return itertools.imap(self.make_val,iter(self._gen()))

  def _gen(self):
    raise ValueError('cannot call this from generic class.  should be overridden')

  def make_val(self,vars):
    raise ValueError('cannot call this from generic class.  should be overridden')

  def _get_block(self):
    """Just read a single block from your current location in _fh"""
    b = self._fh.read(4) # get block size bytes
    #print self._fh.tell()
    if not b: raise StopIteration
    block_size = struct.unpack('<i',b)[0]
    return self._fh.read(block_size)

  def _read_reference_information(self):
    """Reads the reference names and lengths"""
    ref_lengths = {}
    ref_names = []
    for n in range(self._n_ref):
      l_name = struct.unpack('<i',self._fh.read(4))[0]
      name = self._fh.read(l_name).rstrip('\0')
      l_ref = struct.unpack('<i',self._fh.read(4))[0]
      ref_lengths[name] = l_ref
      ref_names.append(name)
    return ref_lengths, ref_names

  def _read_top_header(self):
    """Read the header text and number of reference seqs"""
    magic = self._fh.read(4)
    l_text = struct.unpack('<i',self._fh.read(4))[0]
    header_text = self._fh.read(l_text).rstrip('\0')
    n_ref = struct.unpack('<i',self._fh.read(4))[0]
    return header_text, n_ref

"""BAM File streamer that allows for random access"""
BAMFileOptions = namedtuple('BAMFileOptions',
   ['reference', # reference dictionary
    'blockStart',
    'innerStart',
    'payload'])

class BAMFile(BAMFileGeneric):
  """iterable class to open and access a bam file

  :param filename:
  :param blockStart:
  :param innerStart:
  :param cnt:
  :param reference: dictionary of genomic sequences
  :type filename: string
  :type blockStart: int
  :type innerStart: int
  :type cnt: int
  :type reference: dict()
  """
  def __init__(self,filename,options=None):
    if not options: options = BAMFile.Options()
    super(BAMFile,self).__init__(filename,options)

    self._fh = BGZF(filename)
    self._fetch_headers() #called after ._fh is set

    if self._options.blockStart is not None and self._options.innerStart is not None:
      self._fh.seek(self._options.blockStart,self._options.innerStart)

  @staticmethod
  def Options(**kwargs):
     """ A method for declaring options for the class"""
     construct = BAMFileOptions #IMPORTANT!  Set this
     names = construct._fields
     d = {}
     for name in names: d[name] = None #default values
     for k,v in kwargs.iteritems():
       if k in names: d[k] = v
       else: raise ValueError('Error '+k+' is not a property of these options')
     """Create a set of options based on the inputs"""
     return construct(**d)

  @property
  def reference(self):
    return self._options.reference

  def read_entry(self):
     return self.make_val(self._gen().next())

  def _gen(self):
   while True:
    bstart = self._fh.get_block_start()
    innerstart = self._fh.get_inner_start()
    yield (self._get_block(), self._ref_names,bstart,innerstart)

  def make_val(self,vars):
    data, names, blk, inner = vars
    return BAM(data,names,options = BAM.Options(
               blockStart=blk,
               innerStart=inner,
               reference=self._options.reference,
               header = self.header))
    return vars

  # only get a single
  def fetch_by_coord(self,coord):
    """get a single entry by the coordinate location [blockStart, innerStart]

    .. warning:: creates a new instance of a BAMFile object when maybe the one we had would have worked
    """
    #print coord
    #print self.path
    #b2 = BAMFile(self.path,blockStart=coord[0],innerStart=coord[1],index_obj=self.index,reference=self._reference)
    b2 = BAMFile(self.path,BAMFile.Options(blockStart=coord[0],innerStart=coord[1],reference=self.reference))
    #for bam in b2: print type(bam)
    #print 'hi'
    bam = b2.read_entry()
    b2.close()
    b2 = None
    return bam

  def fetch_starting_at_coord(self,coord):
    #b2 = BAMFile(self.path,blockStart=coord[0],innerStart=coord[1],index_obj=self.index,reference=self._reference)
    """starting at a certain coordinate was supposed to make output

    .. warning:: creates a new instance of a BAMFile object when maybe the one we had would have worked
    """
    b2 = BAMFile(self.path,BAMFile.Options(blockStart=coord[0],innerStart=coord[1],reference=self.reference))
    return b2

"""BAM File streamer that allows for random access"""
BAMFileAltOptions = namedtuple('BAMFileAltOptions',
   ['reference', # reference dictionary
    'payload',
    'mode' # gunzip for subprocess or py for our method seqtools for fastest
   ])

class BAMFileAlt(BAMFileGeneric):
  """iterable class to open and access a bam file

  :param filename:
  :param blockStart:
  :param innerStart:
  :param cnt:
  :param reference: dictionary of genomic sequences
  :type filename: string
  :type blockStart: int
  :type innerStart: int
  :type cnt: int
  :type reference: dict()
  """
  def __init__(self,filename,options=None):
    if not options: options = BAMFileAlt.Options()
    super(BAMFileAlt,self).__init__(filename,options)

    if self._options.mode == 'seqtools':
       self._p0 = open(self._path,'rb')
       self._p1 = Popen('seq-tools bgzf -x -'.split(),stdin=self._p0,stdout=PIPE)
       self._fh = self._p1.stdout
    elif self._options.mode == 'gunzip':
       self._p0 = open(self._path,'rb',10000000)
       self._p1 = Popen('gunzip'.split(),stdin=self._p0,stdout=PIPE,bufsize=10000000)
       self._fh = self._p1.stdout
    elif self._options.mode == 'py':
      self._fh = BGZF(filename)
    else: raise ValueError('Mode '+self._options.mode+" is not valid")

    self._fetch_headers() #called after ._fh is set

  def close(self):
     if self._options.mode == 'gunzip':
       #self._p1.communicate()
       self._p0.close()
     else: self._fh.close()

  @staticmethod
  def Options(**kwargs):
     """ A method for declaring options for the class"""
     construct = BAMFileAltOptions #IMPORTANT!  Set this
     names = construct._fields
     d = {}
     for name in names: d[name] = None #default values
     d['mode'] = 'seqtools'
     for k,v in kwargs.iteritems():
       if k in names: d[k] = v
       else: raise ValueError('Error '+k+' is not a property of these options')
     """Create a set of options based on the inputs"""
     return construct(**d)

  def _gen(self):
   while True:
    yield (self._get_block(), self._ref_names)

  def make_val(self,vars):
    data, names = vars
    return BAM(data,names,options = BAM.Options(reference=self._options.reference))
    return vars


BGZFChunk = namedtuple('BGZFChunk',['block_size','data'])
class BGZF:
  """ Methods adapted from biopython's bgzf.py

  .. warning:: We already have a BGZF class, i winder why we don't put this there

  :param filename:
  :param blockStart:
  :param innerStart:
  :type filename: string
  :type blockStart: int
  :type innerStart: int
  """
  def __init__(self,filename=None,filehandle=None,blockStart=None,innerStart=None,check_crc=False):
    self.path = None
    self._buffer_pos = 0
    self._block_start = 0
    if filename:
       self.path = filename
       self.fh = open(filename,'rb',1000000)
       if blockStart is not None and innerStart is not None: 
         self.seek(blockStart,innerStart)
    elif filehandle:
       self.fh = filehandle
    self.check_crc = check_crc
    self._buffer = self._load_block()
  def close(self):
    self.fh.close()
  def get_block_start(self):
    return self._block_start
  def get_inner_start(self):
    return self._buffer_pos
  def seek(self,blockStart,innerStart):
    self.fh.seek(blockStart)
    self._buffer_pos = 0
    self._buffer = self._load_block()
    self._buffer_pos = innerStart
    """go to this posiiton in the file"""
  def read(self,size):
    """Read this many bytes from where you currently are"""
    done = 0 #number of bytes that have been read so far
    v = ''
    while True:
      if size-done < len(self._buffer.data) - self._buffer_pos:
        v +=  self._buffer.data[self._buffer_pos:self._buffer_pos+(size-done)]
        self._buffer_pos += (size-done)
        #self.pointer += size
        return v
      else: # we need more buffer
        vpart = self._buffer.data[self._buffer_pos:]
        self._buffer = self._load_block()
        v += vpart
        self._buffer_pos = 0
        if len(self._buffer.data)==0: return v
        done += len(vpart)

  def _load_block(self):
    #pointer_start = self.fh.tell()
    if not self.fh: return BGZFChunk(block_size=0,data='')
    self._block_start = self.fh.tell()
    magic = self.fh.read(4)
    #print struct.unpack("<I",magic)[0]
    #print self.get_block_start()
    #print self.get_inner_start()
    if len(magic) < 4:
      #print 'end?'
      #print len(self.fh.read())
      return BGZFChunk(block_size=0,data='')
    gzip_mod_time, gzip_extra_flags, gzip_os,extra_len = struct.unpack("<LBBH",self.fh.read(8))
    pos = 0
    block_size = None
    # need to find the block size
    while pos < extra_len:
      subfield_id = self.fh.read(2)
      subfield_len = struct.unpack("<H",self.fh.read(2))[0]
      subfield_data = self.fh.read(subfield_len)
      pos += subfield_len+4
      if subfield_id == 'BC':
        block_size = struct.unpack("<H",subfield_data)[0]+1
    #block_size is determined
    deflate_size = block_size - 1 - extra_len - 19
    d = zlib.decompressobj(-15)
    data = d.decompress(self.fh.read(deflate_size))+d.flush()
    if not self.check_crc:
      self.fh.read(8)
    else:
      expected_crc = self.fh.read(4)
      expected_size = struct.unpack("<I",self.fh.read(4))[0]
      if expected_size != len(data):
        sys.stderr.write("ERROR unexpected size\n")
        sys.exit()
      crc = zlib.crc32(data)
      if crc < 0:  crc = struct.pack("<i",crc)
      else:  crc = struct.pack("<I",crc)
      if crc != expected_crc:
        sys.stderr.write("ERROR crc fail\n")
        sys.exit()
    return BGZFChunk(block_size=block_size,data=data)
