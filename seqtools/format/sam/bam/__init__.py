"""Classes to work with sam and bam files"""
import struct, zlib, sys, re
from collections import namedtuple
import seqtools.format.sam
from seqtools.format.sam.header import SAMHeader
#from seqtools.sequence import rc
from cStringIO import StringIO
from string import maketrans
_bam_ops = maketrans('012345678','MIDNSHP=X')
_bam_char = maketrans('abcdefghijklmnop','=ACMGRSVTWYHKDBN')
_bam_value_type = {'c':[1,'<b'],'C':[1,'<B'],'s':[2,'<h'],'S':[2,'<H'],'i':[4,'<i'],'I':[4,'<I']}


"""BAM entry (much like a sam line)"""
BAMOptions = namedtuple('BAMOptions',
   ['reference', # reference dictionary
    'header',
    'blockStart',
    'innerStart',
    'payload'])

BAMFields = namedtuple('BAMFields',
   ['qname',
    'flag',
    'rname',
    'pos',
    'mapq',
    'cigar_bytes',
    'rnext',
    'pnext',
    'tlen',
    'seq_bytes',
    'qual_bytes',
    'extra_bytes'])

class BAM(seqtools.format.sam.SAM):
  """Very much like a sam entry but optimized for access from a bam
  Slows down for accessing things that need more decoding like
  sequence, quality, cigar string, and tags

  .. warning:: Having the reference names and the reference sizes we may save some time by not reading the header each time we access the file.  Theres probably a more efficent coarse to go by defining a bamfile object and having a BAMline entry being the extension of sam, and drawing most of this stuff from the bam file

  :param bin_data: byte data for just a single bam entry (seems unnecessary since we have the file)
  :param ref_names: array of refernece names
  :param fileName: the bam file name
  :param blockStart:
  :param innerStart:
  :param ref_lengths: seems unncessary to take the reference lengths because we can always get that from the header
  :param reference:
  :param line_number:
  :type bin_data: bytes
  :type ref_names: list of names
  :type blockStart: where to begin in the file
  :type innerStart: where to begin in the decompressed block
  :type ref_lengths: dict()
  :type reference: dict()
  :type line_number: int
  """
  def __init__(self,bin_data,ref_names,options=None):
    if not options: options = BAM.Options()
    self._options = options
    self._bin_data = bin_data
    self._ref_names = ref_names

    self._auxillary_string = None
    self._bentries = None
    self._line = None
    self._target_range = None
    self._tags = None
    self._seq = None
    self._qual = None
    self._cigar = None
    self._cigar_string = None
    self._alignment_ranges = None #should be accessed by method because of BAM
    return

  #Alignment Ranges are calculated in SAM

  @property
  def bentries(self):
    if self._bentries: return self._bentries 
    self._bentries = _parse_bam_data_block(self._bin_data,self._ref_names) #the binary data and the header is enough to parse 
    return self.bentries

  @staticmethod
  def Options(**kwargs):
     """ A method for declaring options for the class"""
     construct = BAMOptions #IMPORTANT!  Set this
     names = construct._fields
     d = {}
     for name in names: d[name] = None #default values
     for k,v in kwargs.iteritems():
       if k in names: d[k] = v
       else: raise ValueError('Error '+k+' is not a property of these options')
     """Create a set of options based on the inputs"""
     return construct(**d)


  def get_coord(self):
    """get the current coordinate

    :return: [blockStart, innerStart]
    :rtype: list is a pair [int, int]
    """
    return [self._options.blockStart,self._options.innerStart]

  def __str__(self):
    return self.get_sam_line()

  def get_sam_line(self):
    return "\t".join([str(x) for x in
      [self.qname,
       self.flag,
       self.rname,
       self.pos,
       self.mapq,
       self.cigar_string,
       self.rnext,
       self.pnext,
       self.tlen,
       self.seq,
       self.qual,
       self.auxillary_string]])
  ### The getters for all the fields
  @property
  def qname(self): return self.bentries.qname
  @property
  def flag(self): return self.bentries.flag
  @property
  def rname(self): return self.bentries.rname
  @property
  def pos(self): return self.bentries.pos
  @property
  def mapq(self): return self.bentries.mapq

  @property
  def cigar(self): 
    """produce the cigar in list form

    :return: Cigar list of [value (int), type (char)] pairs
    :rtype: list
    """
    if self._cigar: return self._cigar
    v1,v2 = _bin_to_cigar(self.bentries.cigar_bytes)
    self._cigar_string = v2
    self._cigar = v1
    return self._cigar

  @property
  def cigar_string(self):
    if self._cigar_string: return self._cigar_string
    v1,v2 = _bin_to_cigar(self.bentries.cigar_bytes)
    self._cigar_string = v2
    self._cigar = v1
    return self._cigar_string

  @property
  def rnext(self): return self.bentries.rnext  
  @property
  def pnext(self): return self.bentries.pnext  
  @property
  def tlen(self): return self.bentries.tlen
  
  @property
  def seq(self):
    if self._seq: return self._seq
    self._seq = _bin_to_seq(self.bentries.seq_bytes)
    if not self._seq: self._seq = '*'
    return self._seq

  @property
  def qual(self):
    if self._qual: return self._qual
    self._qual = _bin_to_qual(self.bentries.qual_bytes)
    if not self._qual: self._qual = '*'
    return self._qual

  @property
  def tags(self):
    if self._tags: return self._tags
    self._tags, self._auxillary_string = _bin_to_extra(self.bentries.extra_bytes)
    return self._tags

  @property 
  def auxillary_string(self):
    if self._auxillary_string: return self._auxillary_string
    self._tags, self._auxillary_string = _bin_to_extra(self.bentries.extra_bytes)
    return self._auxillary_string

# reference is a dict
"""BAM File streamer"""
BAMFileOptions = namedtuple('BAMFileOptions',
   ['reference', # reference dictionary
    'blockStart',
    'innerStart',
    'payload'])

class BAMFile:
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
    self._options = options
    self._path = filename
    self._fh = BGZF(filename)

    self._line_number = 0 # entry line number ... after header.  starts with 1


    self._header_text, self._n_ref = self._read_top_header()
    self._ref_lengths, self._ref_names = self._read_reference_information()
    self._header = SAMHeader(self._header_text)

    if self._options.blockStart is not None and self._options.innerStart is not None:
      self._fh.seek(blockStart,innerStart)

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
  def header(self):
    return self._header

  def close(self):
    self._fh.close()

  def __iter__(self):
    return self

  def read_entry(self):
    e = self.read_entry2()
    #print e
    #if self._output_range: # check and see if we are past out put range
    #  if not e.is_aligned(): 
    #    e = None
    #  else:
    #    rng2 = e.get_target_range()
    #    if self._output_range.chr != rng2.chr: e = None 
    #    if self._output_range.cmp(rng2) == 1: e = None
    if not e:
      return None
    else: return e

  def next(self):
    e = self.read_entry()
    if not e:
      raise StopIteration
    else: return e

  def read_entry2(self):
    bstart = self._fh.get_block_start()
    innerstart = self._fh.get_inner_start()
    b = self._fh.read(4) # get block size bytes
    if not b: return None
    block_size = struct.unpack('<i',b)[0]
    #print 'block_size '+str(block_size)
    self._line_number += 1
    bam = BAM(self._fh.read(block_size),self._ref_names,options = BAM.Options(
      blockStart=bstart,innerStart=innerstart,header=self.header))
    return bam

  #def _set_output_range(self,rng):
  #  self._output_range = rng
  #  return

  # only get a single
  def fetch_by_coord(self,coord):
    """get a single entry by the coordinate location [blockStart, innerStart]

    .. warning:: creates a new instance of a BAMFile object when maybe the one we had would have worked
    """
    #b2 = BAMFile(self.path,blockStart=coord[0],innerStart=coord[1],index_obj=self.index,reference=self._reference)
    b2 = BAMFile(self.path,blockStart=coord[0],innerStart=coord[1],reference=self._reference)
    bam = b2.read_entry()
    b2.close()
    b2 = None
    return bam

  def fetch_starting_at_coord(self,coord):
    #b2 = BAMFile(self.path,blockStart=coord[0],innerStart=coord[1],index_obj=self.index,reference=self._reference)
    """starting at a certain coordinate was supposed to make output

    .. warning:: creates a new instance of a BAMFile object when maybe the one we had would have worked
    """
    b2 = BAMFile(self.path,blockStart=coord[0],innerStart=coord[1],reference=self._reference)
    return b2

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

def _parse_bam_data_block(bin_in,ref_names):
  data = StringIO(bin_in)
  rname_num = struct.unpack('<i',data.read(4))[0]
  v_rname = ref_names[rname_num] #refID to check in ref names
  v_pos = struct.unpack('<i',data.read(4))[0] + 1 #POS
  bin_mq_nl = struct.unpack('<I',data.read(4))[0]
  bin =  bin_mq_nl >> 16 
  v_mapq = (bin_mq_nl & 0xFF00) >> 8 #mapq
  l_read_name = bin_mq_nl & 0xFF #length of qname
  flag_nc = struct.unpack('<I',data.read(4))[0] #flag and n_cigar_op
  v_flag = flag_nc >> 16
  n_cigar_op = flag_nc & 0xFFFF
  l_seq = struct.unpack('<i',data.read(4))[0]
  rnext_num = struct.unpack('<i',data.read(4))[0]
  if rnext_num == -1:
    v_rnext = '*'
  else:
    v_rnext = ref_names[rnext_num] #next_refID in ref_names
  v_pnext = struct.unpack('<i',data.read(4))[0]+1 #pnext
  tlen = struct.unpack('<i',data.read(4))[0]
  v_tlen = tlen
  v_qname = data.read(l_read_name).rstrip('\0') #read_name or qname
  #print 'n_cigar_op '+str(n_cigar_op)
  v_cigar_bytes = data.read(n_cigar_op*4)
  #print 'cigar bytes '+str(len(v['cigar_bytes']))
  v_seq_bytes = data.read((l_seq+1)/2)
  v_qual_bytes = data.read(l_seq)
  v_extra_bytes = data.read()
  #last second tweak
  if v_rnext == v_rname: v_rnext = '='
  return BAMFields(
     v_qname,
     v_flag,
     v_rname,
     v_pos,
     v_mapq,
     v_cigar_bytes,
     v_rnext,
     v_pnext,
     v_tlen,
     v_seq_bytes,
     v_qual_bytes,
     v_extra_bytes)

def _bin_to_qual(qual_bytes):
  if len(qual_bytes) == 0: return '*'
  if struct.unpack('<B',qual_bytes[0])[0] == 0xFF: return '*'
  #print qual_bytes
  #try:
  qual = ''.join([chr(struct.unpack('<B',x)[0]+33) for x in qual_bytes])
  #except:
  #  return '*'
  return qual

def _bin_to_seq(seq_bytes):
  if len(seq_bytes) == 0: return None
  global _bam_char
  #print len(seq_bytes)
  seq = ''.join([''.join([''.join([chr(z+97).translate(_bam_char) for z in  [y>>4,y&0xF]]) for y in struct.unpack('<B',x)]) for x in seq_bytes]).rstrip('=')
  return seq

def _bin_to_cigar(cigar_bytes):
  global _bam_ops
  if len(cigar_bytes) == 0: return [[],'*']
  cigar_packed = [struct.unpack('<I',x)[0] for x in \
             [cigar_bytes[i:i+4] for i in range(0,len(cigar_bytes),4)]]
  cigar_array = [[c >> 4, str(c &0xF).translate(_bam_ops)] for c in cigar_packed]
  cigar_seq = ''.join([''.join([str(x[0]),x[1]]) for x in cigar_array])
  return [cigar_array,cigar_seq]

#Pre all the reamining bytes of an entry
#Post an array of 
# 1. A dict keyed by Tag with {'type':,'value':} where value is a string unless type is i
# 2. A string of the remainder
def _bin_to_extra(extra_bytes):
  global _bam_value_type
  extra = StringIO(extra_bytes)
  tags = {}
  rem = ''
  while extra.tell() < len(extra_bytes):
    tag = extra.read(2)
    val_type = extra.read(1)
    if val_type == 'Z':
      rem += tag+':'
      rem += val_type+':'
      p = re.compile('([!-~])')
      m = p.match(extra.read(1))
      vre = ''
      while m:
        vre += m.group(1)
        c = extra.read(1)
        #print c
        m = p.match(c)
      rem += vre+"\t"
      tags[tag] = {'type':val_type,'value':vre}
    elif val_type == 'A':
      rem += tag+':'
      rem += val_type+':'
      vre = extra.read(1)
      rem += vre+"\t"      
      tags[tag] = {'type':val_type,'value':vre}      
    elif val_type in _bam_value_type:
      rem += tag+':'
      rem += 'i'+':'
      val = struct.unpack(_bam_value_type[val_type][1],extra.read(_bam_value_type[val_type][0]))[0]
      rem += str(val)+"\t"
      tags[tag] = {'type':val_type,'value':val}
    elif val_type == 'B':
      sys.sterr.write("WARNING array not implmented\n")
      continue
      rem += tag+':'
      rem += val_type+':'
      array_type = _bam_value_type[extra.read(1)]
      element_count = struct.unpack('<I',extra.read(4))[0]
      array_bytes = extra.read(element_count*_bam_value_type[array_type][0])
      for by in [array_bytes[i:i+_bam_value_type[array_type][0]] for i in range(0,len(array_bytes),_bam_value_type[array_type][0])]:
        aval = struct.unpack(_bam_value_type[array_type][1],by)
  return [tags,rem.rstrip("\t")]

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
  def __init__(self,filename,blockStart=None,innerStart=None):
    self.path = filename
    self.fh = open(filename,'rb')
    if blockStart: self.fh.seek(blockStart)
    self._block_start = 0
    self._buffer = self._load_block()
    self._buffer_pos = 0
    if innerStart: self._buffer_pos = innerStart
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
    if len(magic) < 4:
      #print 'end?'
      #print len(self.fh.read())
      return BGZFChunk(block_size=0,data='')
    gzip_mod_time, gzip_extra_flags, gzip_os,extra_len = struct.unpack("<LBBH",self.fh.read(8))
    pos = 0
    block_size = None
    #get block_size
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
