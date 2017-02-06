"""Classes to work with sam and bam files"""
import struct, zlib, sys, re, os, gzip, random
import seqtools.align
import seqtools.format.sam
from seqtools.sequence import rc
from cStringIO import StringIO
from string import maketrans
from seqtools.range import GenomicRange
from subprocess import Popen, PIPE
_bam_ops = maketrans('012345678','MIDNSHP=X')
_bam_char = maketrans('abcdefghijklmnop','=ACMGRSVTWYHKDBN')
_bam_value_type = {'c':[1,'<b'],'C':[1,'<B'],'s':[2,'<h'],'S':[2,'<H'],'i':[4,'<i'],'I':[4,'<I']}
_sam_cigar_target_add = re.compile('[M=XDN]$')

class BAM(samtools.format.sam.SAM):
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
  def __init__(self,bin_data,ref_names,fileName=None,blockStart=None,innerStart=None,ref_lengths=None,reference=None,line_number=None):
    part_dict = _parse_bam_data_block(bin_data,ref_names)
    #self._bamfileobj = bamfileobj #this is most like our parent
    self._line = None
    self._line_number = line_number # the line number in the bam file
    self._reference = reference
    self._target_range = None
    self._alignment_ranges = None #should be accessed by method because of BAM
    self._ref_lengths = ref_lengths
    self._file_position = {'fileName':fileName,'blockStart':blockStart,'innerStart':innerStart} # The most special information about the bam
    self._private_values = BAM.PrivateValues() # keep from accidently accessing some variables other than by methods
    self._private_values.set_entries_dict(part_dict)
    #self._set_alignment_ranges()
    return

  def get_alignment_ranges(self):
    """return the basics for defining an alignment"""
    if not self._alignment_ranges:
      self._set_alignment_ranges()
    return self._alignment_ranges

  def get_line_number(self):
    return self._line_number
  def get_target_length(self):
    """length of the entire chromosome"""
    return self._ref_lengths[self.value('rname')]
  def get_filename(self):
    return self._file_position['fileName']
  def get_coord(self):
    """get the current coordinate

    :return: [blockStart, innerStart]
    :rtype: list is a pair [int, int]
    """
    return [self._file_position['blockStart'],self._file_position['innerStart']]
  def get_block_start(self):
    return self._file_position['blockStart']
  def get_inner_start(self):
    return self._file_position['innerStart']
  def get_file_position_string(self):
    return 'fileName: '+self._file_position['fileName']+" "\
           'blockStart: '+str(self._file_position['blockStart'])+" "\
           'innerStart: '+str(self._file_position['innerStart'])
  def get_tag(self,key):
    """retrieve the value of a single tag by its key.  

    .. warning:: Not sure if it accommodates multiple of the same keys"""
    cur = self._private_values.get_tags()
    if not cur:
      v1,v2 = _bin_to_extra(self.value('extra_bytes'))
      self._private_values.set_tags(v1) #keep the cigar array in a special palce
      self._private_values.set_entry('remainder',v2)
    return self._private_values.get_tags()[key]['value']
  def get_cigar(self): 
    """produce the cigar in list form

    :return: Cigar list of [value (int), type (char)] pairs
    :rtype: list
    """
    cur = self._private_values.get_cigar()
    if not cur:
      v1,v2 = _bin_to_cigar(self.value('cigar_bytes'))
      self._private_values.set_cigar(v1) #keep the cigar array in a special palce
      self._private_values.set_entry('cigar',v2)
    return self._private_values.get_cigar()

  def value(self,key):
    """Access basic attributes of BAM by key"""
    if not self._private_values.is_entry_key(key):
      if key == 'seq':
        v = _bin_to_seq(self.value('seq_bytes'))
        if not v: v = '*'
        self._private_values.set_entry('seq',v)
        return v
      elif key == 'qual':
        v = _bin_to_qual(self.value('qual_bytes'))
        if not v: v = '*'
        self._private_values.set_entry('qual',v)
        return v
      elif key == 'cigar':
        v1,v2 = _bin_to_cigar(self.value('cigar_bytes'))
        self._private_values.set_cigar(v1) #keep the cigar array in a special palce
        self._private_values.set_entry('cigar',v2)
        return v2
      elif key == 'remainder':
        v1,v2 = _bin_to_extra(self.value('extra_bytes'))
        self._private_values.set_tags(v1) #keep the cigar array in a special palce
        self._private_values.set_entry('remainder',v2)
        return v2
    return self._private_values.get_entry(key)


# reference is a dict
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
  def __init__(self,filename,blockStart=None,innerStart=None,cnt=None,reference=None):
    self.path = filename
    self._reference = reference # dict style accessable reference
    self.fh = BGZF(filename)
    self._line_number = 0 # entry line number ... after header.  starts with 1
    # start reading the bam file
    self.header_text = None
    self._header = None
    self.n_ref = None
    self._read_top_header()
    self.ref_names = []
    self.ref_lengths = {}
    self._output_range = None
    #self.index = index_obj
    self._read_reference_information()
    # prepare for specific work
    if self.path and blockStart is not None and innerStart is not None:
      self.fh.seek(blockStart,innerStart)
      #if self.index:
      #  lnum = self.index.get_coord_line_number([blockStart,innerStart])
      #  if lnum:
      #    self._line_number = lnum-1 #make it zero indexed
  def close(self):
    self.fh.close()
  # return a string that is the header
  def get_header(self):
    if not self._header:
      self._header = SAMHeader(self.header_text)
      return self._header
    return self._header

  def __iter__(self):
    return self

  def read_entry(self):
    e = self.read_entry2()
    #print e
    if self._output_range: # check and see if we are past out put range
      if not e.is_aligned(): 
        e = None
      else:
        rng2 = e.get_target_range()
        if self._output_range.chr != rng2.chr: e = None 
        if self._output_range.cmp(rng2) == 1: e = None
    if not e:
      return None
    else: return e

  def next(self):
    e = self.read_entry()
    if not e:
      raise StopIteration
    else: return e

  def read_entry2(self):
    bstart = self.fh.get_block_start()
    innerstart = self.fh.get_inner_start()
    b = self.fh.read(4) # get block size bytes
    if not b: return None
    block_size = struct.unpack('<i',b)[0]
    #print 'block_size '+str(block_size)
    self._line_number += 1
    bam = BAM(self.fh.read(block_size),self.ref_names,fileName=self.path,blockStart=bstart,innerStart=innerstart,ref_lengths=self.ref_lengths,reference=self._reference,line_number = self._line_number)
    return bam

  def _set_output_range(self,rng):
    self._output_range = rng
    return

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
    for n in range(self.n_ref):
      l_name = struct.unpack('<i',self.fh.read(4))[0]
      name = self.fh.read(l_name).rstrip('\0')
      l_ref = struct.unpack('<i',self.fh.read(4))[0]
      self.ref_lengths[name] = l_ref
      self.ref_names.append(name)
  def _read_top_header(self):
    magic = self.fh.read(4)
    l_text = struct.unpack('<i',self.fh.read(4))[0]
    self.header_text = self.fh.read(l_text).rstrip('\0')
    self.n_ref = struct.unpack('<i',self.fh.read(4))[0]

def _parse_bam_data_block(bin_in,ref_names):
  v = {}
  data = StringIO(bin_in)
  rname_num = struct.unpack('<i',data.read(4))[0]
  v['rname'] = ref_names[rname_num] #refID to check in ref names
  v['pos'] = struct.unpack('<i',data.read(4))[0] + 1 #POS
  bin_mq_nl = struct.unpack('<I',data.read(4))[0]
  bin =  bin_mq_nl >> 16 
  v['mapq'] = (bin_mq_nl & 0xFF00) >> 8 #mapq
  l_read_name = bin_mq_nl & 0xFF #length of qname
  flag_nc = struct.unpack('<I',data.read(4))[0] #flag and n_cigar_op
  v['flag'] = flag_nc >> 16
  n_cigar_op = flag_nc & 0xFFFF
  l_seq = struct.unpack('<i',data.read(4))[0]
  rnext_num = struct.unpack('<i',data.read(4))[0]
  if rnext_num == -1:
    v['rnext'] = '*'
  else:
    v['rnext'] = ref_names[rnext_num] #next_refID in ref_names
  v['pnext'] = struct.unpack('<i',data.read(4))[0]+1 #pnext
  tlen = struct.unpack('<i',data.read(4))[0]
  v['tlen'] = tlen
  v['qname'] = data.read(l_read_name).rstrip('\0') #read_name or qname
  #print 'n_cigar_op '+str(n_cigar_op)
  v['cigar_bytes'] = data.read(n_cigar_op*4)
  #print 'cigar bytes '+str(len(v['cigar_bytes']))
  v['seq_bytes'] = data.read((l_seq+1)/2)
  v['qual_bytes'] = data.read(l_seq)
  v['extra_bytes'] = data.read()
  #last second tweak
  if v['rnext'] == v['rname']: v['rnext'] = '='
  return v

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
    #self.pointer = 0
    #holds block_size and data
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
      if size-done < len(self._buffer['data']) - self._buffer_pos:
        v +=  self._buffer['data'][self._buffer_pos:self._buffer_pos+(size-done)]
        self._buffer_pos += (size-done)
        #self.pointer += size
        return v
      else: # we need more buffer
        vpart = self._buffer['data'][self._buffer_pos:]
        self._buffer = self._load_block()
        v += vpart
        self._buffer_pos = 0
        if len(self._buffer['data'])==0: return v
        done += len(vpart)

  def _load_block(self):
    #pointer_start = self.fh.tell()
    if not self.fh: return {'block_size':0,'data':''}
    self._block_start = self.fh.tell()
    magic = self.fh.read(4)
    if len(magic) < 4:
      #print 'end?'
      #print len(self.fh.read())
      return {'block_size':0,'data':''}
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
    return {'block_size':block_size, 'data':data}
