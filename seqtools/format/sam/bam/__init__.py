"""Classes to work with sam and bam files"""
import struct, zlib, sys, re, itertools
from collections import namedtuple
import seqtools.format.sam
from seqtools.format.sam.header import SAMHeader
#from seqtools.sequence import rc
#from cStringIO import StringIO
from io import StringIO, BytesIO
#from string import maketrans
_bam_ops = str.maketrans('012345678','MIDNSHP=X')
_bam_char = str.maketrans('abcdefghijklmnop','=ACMGRSVTWYHKDBN')
_bam_value_type = {'c':[1,'<b'],'C':[1,'<B'],'s':[2,'<h'],'S':[2,'<H'],'i':[4,'<i'],'I':[4,'<I']}

from seqtools.format.sam import TagDatum, CIGARDatum, check_flag

class BAMEntries:
   """BAM entry data in same format as SAM entries tuple"""
   def __init__(self,binary_data):
      self._bentries = binary_data
      self._cigar_string = None
      self._cigar_array = None
      self._seq = None
      self._qual = None
      self._optional_fields = None
   def is_aligned(self):
      return not check_flag(self.flag,0x4)
   ### The getters for all the fields
   @property
   def qname(self): return self._bentries.qname
   @property
   def flag(self): return self._bentries.flag
   @property
   def rname(self): 
      if not self.is_aligned(): return '*'
      return self._bentries.rname
   @property
   def pos(self):
      if self._bentries.pos == 0: return None 
      return self._bentries.pos
   @property
   def mapq(self): 
      if self._bentries.mapq == 255: return None
      return self._bentries.mapq

   @property
   def cigar(self):
      if self._cigar_string:
         if self._cigar_string == '*': return None 
         return self._cigar_string
      v1,v2 = _bin_to_cigar(self._bentries.cigar_bytes)
      if not v2: self._cigar_string = '*'
      else: self._cigar_string = v2
      self._cigar_array = v1
      return self.cigar

   def get_cigar_array(self):
      c = self.cigar
      return self._cigar_array

   @property
   def rnext(self): return self._bentries.rnext  
   @property
   def pnext(self): return self._bentries.pnext  
   @property
   def tlen(self): return self._bentries.tlen
  
   @property
   def seq(self):
    if self._seq:
       if self._seq == '*': return None
       return self._seq
    self._seq = _bin_to_seq(self._bentries.seq_bytes)
    if not self._seq: self._seq = '*'
    return self.seq
    #if self._seq == '*': return None
    #return self._seq

   @property
   def qual(self):
    if self._qual:
       if self._qual == '*': return None 
       return self._qual
    self._qual = _bin_to_qual(self._bentries.qual_bytes)
    if not self._qual: self._qual = '*'
    return self.qual
    #if self._qual == '*': return None 
    #return self._qual

   @property 
   def optional_fields(self):
      if self._optional_fields: 
         if self._optional_fields == '': return None
         return self._optional_fields
      self._tags, self._optional_fields = _bin_to_extra(self._bentries.extra_bytes)
      if not self._optional_fields: self._optional_fields = ''
      return self.optional_fields

   def get_tags(self):
      o = self.optional_fields
      return self._tags

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
    self.bentries = _parse_bam_data_block(self._bin_data,self._ref_names) #the binary data and the header is enough to parse 
    self._entries = BAMEntries(self.bentries)
    self._line = None
    self._target_range = None
    self._alignment_ranges = None #should be accessed by method because of BAM
    return

  @property
  def entries(self): return self._entries

  #Alignment Ranges are calculated in SAM

  @staticmethod
  def Options(**kwargs):
     """ A method for declaring options for the class"""
     construct = BAMOptions #IMPORTANT!  Set this
     names = construct._fields
     d = {}
     for name in names: d[name] = None #default values
     for k,v in kwargs.items():
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

  @property
  def blockStart(self):
    """Maybe none if not set"""
    return self._options.blockStart
  @property
  def innerStart(self):
    """Maybe none if not set"""
    return self._options.innerStart

  def __str__(self):
    return self.sam_line

  @property
  def sam_line(self):
     out = self.entries.qname + "\t"
     out += str(self.entries.flag) + "\t"
     if self.entries.rname is None:
        out += '*' + "\t"
     else:
        out += self.entries.rname + "\t" 
     if self.entries.pos is None:
        out += '0' + "\t"
     else:
        out += str(self.entries.pos) + "\t"
     if self.entries.mapq is None:
        out += '255'+"\t"
     else:
        out += str(self.entries.mapq) + "\t"

     if self.entries.cigar is None:
        out += '*'+"\t"
     else:
        out += self.entries.cigar+"\t"

     if self.entries.rnext is None:
        out += '*'+"\t"
     else:
        out += self.entries.rnext+"\t"

     if self.entries.pnext is None:
        out += '0'+"\t"
     else:
        out += str(self.entries.pnext)+"\t"
     if self.entries.tlen is None:
        out += '0'+"\t"
     else:
        out += str(self.entries.tlen)+"\t"
     if self.entries.seq is None:
        out += '*'+"\t"
     else:
        out += self.entries.seq+"\t"
     if self.entries.qual is None:
        out += '*'
     else:
        out += self.entries.qual
     if self.entries.optional_fields:
        out += "\t"+self.entries.optional_fields
     return out

  @property
  def cigar_array(self): 
    """produce the cigar in list form

    :return: Cigar list of [value (int), type (char)] pairs
    :rtype: list
    """
    return self.entries.get_cigar_array()

  @property
  def tags(self):
    return self.entries.get_tags()

def _parse_bam_data_block(bin_in,ref_names):
  data = BytesIO(bin_in)
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
  #global _bam_char
  #print len(seq_bytes)
  seq = ''.join([''.join(
           [''.join([chr(z+97).translate(_bam_char) for z in  [y>>4,y&0xF]]) for y in struct.unpack('<B',x)]) for x in seq_bytes]).rstrip('=')
  return seq

def _bin_to_cigar(cigar_bytes):
  #global _bam_ops
  if len(cigar_bytes) == 0: return [[],'*']
  cigar_packed = [struct.unpack('<I',x)[0] for x in \
             [cigar_bytes[i:i+4] for i in range(0,len(cigar_bytes),4)]]
  cigar_array = [CIGARDatum(c >> 4, str(c &0xF).translate(_bam_ops)) for c in cigar_packed]
  cigar_seq = ''.join([''.join([str(x[0]),x[1]]) for x in cigar_array])
  return [cigar_array,cigar_seq]

def _bin_to_extra(extra_bytes):
  """Pre all the reamining bytes of an entry
  Post an array of 
   1. A dict keyed by Tag with {'type':,'value':} where value is a string unless type is i
   2. A string of the remainder
  """
  #global _bam_value_type
  extra = BytesIO(extra_bytes)
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
      tags[tag] = TagDatum(val_type, vre)
      #tags[tag] = {'type':val_type,'value':vre}
    elif val_type == 'A':
      rem += tag+':'
      rem += val_type+':'
      vre = extra.read(1)
      rem += vre+"\t"      
      tags[tag] = TagDatum(val_type, vre)
      #tags[tag] = {'type':val_type,'value':vre}      
    elif val_type in _bam_value_type:
      rem += tag+':'
      rem += 'i'+':'
      val = struct.unpack(_bam_value_type[val_type][1],extra.read(_bam_value_type[val_type][0]))[0]
      rem += str(val)+"\t"
      tags[tag] = TagDatum(val_type, val)
      #tags[tag] = {'type':val_type,'value':val}
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
    self.fh = open(filename,'rb',1000000)
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
