"""Classes to work with sam and files"""
import sys, re, os, itertools
from collections import namedtuple

import seqtools.align
from seqtools.sequence import rc
from seqtools.range import GenomicRange
from seqtools.format.sam.header import SAMHeader
import seqtools.stream
from tempfile import gettempdir
from subprocess import Popen, PIPE
_sam_cigar_target_add = re.compile('[M=XDN]$')

TagDatum = namedtuple('TAGDatum',
   ['type','value'])
CIGARDatum = namedtuple('CIGARDatum',
   ['value','op'])

"""SAM line options that are not absolutely necessary for a sam line"""
SAMOptions = namedtuple('SAMOptions',
   ['reference', # reference dictionary
    'reference_lengths', # lengths of chromosomes in reference dictionary
    'header',
    'payload'])

SAMFields = namedtuple('SAMFields',
   ['qname',
    'flag',
    'rname',
    'pos',
    'mapq',
    'cigar',
    'rnext',
    'pnext',
    'tlen',
    'seq',
    'qual',
    'optional_fields'])

class SAM(seqtools.align.Alignment):
  """Class to define the SAM format line (not header).
  :param line:  a SAM line
  :param options:
  :param options.sam_header: an object representing the sam header
  :param options.reference:  a reference dict of sequences
  :param options.reference_lengths: reference dict of sequence lengths
  :type line: string
  :type options: namedtuple
  :type options.sam_header: SAMHeader object
  :type options.reference: dict() 
  :type options.reference_lengths: dict() - dictionary of chromosome keyed lengths
  """
  def __init__(self,line,options=None):
    if not options: options = SAM.Options()
    self._options = options
    self._line = line.rstrip()
    self._cigar = None # Cache for the cigar
    self._tags = None # Cache fro tags
    self._target_range = None
    self._alignment_ranges = None
    self._entries = None
    if self.is_aligned():
       super(SAM,self).__init__(options)
    return

  @property
  def sam_line(self):
     return self._line

  @property
  def entries(self):
     if self._entries: return self._entries
     self._entries = self._parse_sam_line()
     return self._entries
  @property
  def alignment_ranges(self):
     """Put the heavy alignment ranges calculation called on demand and then cached"""
     if not self.is_aligned(): raise ValueError("you can't get alignment ranges from something that didn't align")
     if self._alignment_ranges: return self._alignment_ranges
     self._alignment_ranges = self._get_alignment_ranges()
     return self._alignment_ranges

  @staticmethod
  def Options(**kwargs):
     """ A method for declaring options for the class"""
     construct = SAMOptions #IMPORTANT!  Set this
     names = construct._fields
     d = {}
     for name in names: d[name] = None #default values
     for k,v in kwargs.iteritems():
       if k in names: d[k] = v
       else: raise ValueError('Error '+k+' is not a property of these options')
     """Create a set of options based on the inputs"""
     return construct(**d)

  def _parse_sam_line(self):
    f = self._line.rstrip().split("\t")
    qname = f[0]
    flag = int(f[1])
    rname = f[2]
    if rname == '*': rname = None
    pos = int(f[3])
    if pos == 0: pos = None
    mapq = int(f[4])
    if mapq == 255: mapq = None 
    cigar = f[5]
    if f[5] == '*': cigar = None
    rnext = f[6]
    if f[6] == '*': rnext = None
    pnext = int(f[7])
    if pnext == 0: pnext = None 
    tlen = int(f[8])
    if tlen == 0: tlen = None 
    seq = None
    if f[9] != '*': seq = f[9]
    qual = None
    if f[10] != '*': qual = f[10]
    optional_fields = None
    if len(f) > 11:
       optional_fields = "\t".join(f[11:])

    return  SAMFields(
       qname, #qname
       flag, #flag
       rname, #rname
       pos, #pos
       mapq, #mapq
       cigar, #cigar
       rnext, #rnext
       pnext, #pnext
       tlen, #tlen
       seq, #seq
       qual, #qual
       optional_fields)

  def __str__(self):
    """The only way to set a SAM line is with a SAM line so we have this"""
    return self._line

  def get_aligned_bases_count(self):
     if not self.is_aligned(): return 0
     return super(SAM,self).get_aligned_bases_count()

  @property
  def target_sequence_length(self):
    """ Get the length of the target sequence.  length of the entire chromosome

    throws an error if there is no information available

    :return: length
    :rtype: int
    """
    if not self.is_aligned():
      raise ValueError("no length for reference when read is not not aligned")
    if self.entries.tlen: return self.entries.tlen #simplest is if tlen is set
    if self.header:
      if self.entries.rname in self.header.sequence_lengths:
        return self.header.sequence_lengths[self.entries.rname]
    elif self.reference:
      return len(self.reference[self.entries.rname])
    else:
      raise ValueError("some reference needs to be set to go from psl to bam\n")
    raise ValueError("No reference available")

  @property
  def query_sequence(self):
    """ Overrides align. corrects orientation with reverse complement if on negative strand

    .. warning:: this returns the full query sequence, not just the aligned portion, but i also does not include hard clipped portions (only soft clipped)
    """
    if not self.entries.seq: return None
    if self.check_flag(0x10): return rc(self.entries.seq)
    return self.entries.seq

  @property
  def query_quality(self):
    """ Overrides align 

    .. warning:: this returns the full query quality, not just the aligned portion
    """
    if not self.entries.qual: return None
    if self.entries.qual == '*': return None
    if self.check_flag(0x10): return self.entries.qual[::-1]
    return self.entries.qual

  @property
  def query_sequence_length(self):
    """ does not include hard clipped"""
    if self.entries.seq: return len(self.entries.seq)
    if not self.entries.cigar:
       raise ValueError('Cannot give a query length if no cigar and no query sequence are present')
    return sum([x[0] for x in self.cigar_array if re.match('[MIS=X]',x[1])])

  @property
  def original_query_sequence_length(self):
    """Similar to get_get_query_sequence_length, but it also includes
    hard clipped bases
    if there is no cigar, then default to trying the sequence

    :return: the length of the query before any clipping
    :rtype: int
    """
    if not self.is_aligned() or not self.entries.cigar:
      return self.query_sequence_length # take the naive approach 
    # we are here with something aligned so take more intelligent cigar apporach  
    return sum([x[0] for x in self.cigar_array if re.match('[HMIS=X]',x[1])])

  @property
  def actual_original_query_range(self):
    """ This accounts for hard clipped bases 
    and a query sequence that hasnt been reverse complemented

    :return: the range covered on the original query sequence
    :rtype: GenomicRange
    """
    l = self.original_query_sequence_length
    a = self.alignment_ranges
    qname = a[0][1].chr
    qstart = a[0][1].start
    qend = a[-1][1].end
    #rng = self.get_query_range()
    start = qstart
    end = qend
    if self.strand == '-':
      end = l-(qstart-1)
      start = 1+l-(qend)
    return GenomicRange(qname,start,end,dir=self.strand)

  @property
  def strand(self):
    """ Overrides parent to get direction from the flag

    :return: strand/direction + or -
    :rtype: char
    """
    if self.check_flag(0x10): return '-'
    return '+'

  def get_SAM(self):
    """Override parent to just return itself"""
    return self

  def _get_alignment_ranges(self):
    """A key method to extract the alignment data from the line"""
    if not self.is_aligned(): return None
    alignment_ranges = []
    cig = [x[:] for x in self.cigar_array]
    target_pos = self.entries.pos
    query_pos = 1
    while len(cig) > 0:
      c = cig.pop(0)
      if re.match('[S]$',c[1]): # hard or soft clipping
        query_pos += c[0]
      elif re.match('[ND]$',c[1]): # deleted from reference
        target_pos += c[0]
      elif re.match('[I]$',c[1]): # insertion to the reference
        query_pos += c[0]
      elif re.match('[MI=X]$',c[1]): # keep it
        t_start = target_pos
        q_start = query_pos
        target_pos += c[0]
        query_pos += c[0]
        t_end = target_pos-1
        q_end = query_pos-1
        alignment_ranges.append([GenomicRange(self.entries.rname,t_start,t_end),GenomicRange(self.entries.qname,q_start,q_end)])
    return alignment_ranges

  @property
  def range(self):
    """ Necessary property for doing a locus stream
    For the context of a SAM file we set this to be the target range

    :return: target range
    :rtype: GenomicRange
    """
    return self.target_range

  @property
  def target_range(self):
    """Get the range on the target strand

    :return: target range
    :rtype: GenomicRange
    """
    if not self.is_aligned(): return None
    if self._target_range: return self._target_range # check cache
    global _sam_cigar_target_add
    tlen = sum([x[0] for x in self.cigar_array if _sam_cigar_target_add.match(x[1])])
    self._target_range = GenomicRange(self.entries.rname,self.entries.pos,self.entries.pos+tlen-1)
    return self._target_range

  def check_flag(self,inbit):
    if self.entries.flag & inbit: return True
    return False

  def is_aligned(self):
    """Basic but very important since we are an alignment class and not everything is aligned"""
    return not self.check_flag(0x4)

  def get_line(self):
    """assemble the line if its not there yet

    .. warning:: this should probably not exist if the constructor takes a line
    """
    return self._line
  @property
  def header(self):
    return self._options.header

  ##### PUT GETTERS FOR MAIN PROPERTIES HERE ####
  #@property
  #def qname(self): return self.entries.qname
  #@property
  #def flag(self): return self.entries.flag
  #@property
  #def rname(self): return self.entries.rname
  #@property
  #def pos(self): return self.entries.pos
  #@property
  #def mapq(self): return self.entries.mapq

  @property
  def cigar_array(self):
     """cache this one to speed things up a bit"""
     if self._cigar: return self._cigar
     self._cigar = [CIGARDatum(int(m[0]),m[1]) for m in re.findall('([0-9]+)([MIDNSHP=X]+)',self.entries.cigar)]
     return self._cigar

  #@property
  #def cigar_string(self):
  #  """don't get at this one from the entries because it differes in bam"""
  #  return self.entries.cigar

  #@property
  #def rnext(self): return self.entries.rnext
  #@property
  #def pnext(self): return self.entries.pnext
  #@property
  #def tlen(self): return self.entries.tlen

  #@property
  #def seq(self):
  #  """Need to use this to access seq because of BAM"""
  #  return self.entries.seq

  #@property
  #def qual(self):
  #  """Need to use this to access seq because of BAM"""
  #  return self.entries.qual

  @property
  def tags(self):
     """Access the auxillary data here"""
     if self._tags: return self._tags
     tags = {}
     if not tags: return {}
     for m in [[y.group(1),y.group(2),y.group(3)] for y in [re.match('([^:]{2,2}):([^:]):(.+)$',x) for x in self.entries.optional_fields.split("\t")]]:
        if m[1] == 'i': m[2] = int(m[2])
        elif m[1] == 'f': m[2] = float(m[2])
        tags[m[0]] = TAGDatum(m[1],m[2])
     self._tags = tags
     return self._tags

  #@property
  #def auxillary_string(self):
  #  v = self._line.rstrip().split("\t")
  #  if len(v) > 11: return v[11:]
  #  return None

SAMGeneratorOptions = namedtuple('SAMGeneratorOptions',
   ['buffer_size',
    'reference'])

class SAMGenerator(seqtools.stream.BufferedLineGenerator):
   """Generic class for a streaming SAM data

      Borrows the _gen fucntion from BufferedLineGenerator to shorten things a bit
   """
   def __init__(self,stream,options=None):
      if not options: options = SAMGenerator.Options()
      self._options = options
      self._stream = stream
      self._buffer = ''
      self._header_text = ''
      # start falling through the header
      done_header = False
      while True:
         data = self._stream.read(self._options.buffer_size)
         if not data: break
         self._buffer += data
         last = 0
         for m in re.finditer('\n',data):
            ln = data[last:m.start()+1]
            if is_header(ln): 
               self._header_text += ln
            else:
               done_header = True
            last = m.start()+1
         if done_header: break
      #leave the unprocessed portion of the buffer
      self._buffer = self._buffer[len(self._header_text):]
      #After this point it is a normal buffered line generator, using
      #the SAM mapping function

   @staticmethod
   def Options(**kwargs):
      """ A method for declaring options for the class"""
      construct = SAMGeneratorOptions #IMPORTANT!  Set this
      names = construct._fields
      d = {}
      for name in names: d[name] = None #default values
      """set defaults here"""
      d['buffer_size'] = 10000000
      for k,v in kwargs.iteritems():
         if k in names: d[k] = v
         else: raise ValueError('Error '+k+' is not a property of these options')
      """Create a set of options based on the inputs"""
      return construct(**d)

   def __iter__(self):
      return itertools.imap(self.make_sam,self._gen())

   def make_sam(self,line):
      return SAM(line,SAM.Options(reference=self._options.reference,
                                  header=self.header))

   def has_header(self):
      if len(self._header_text) > 0: return True
      return False

   @property
   def header(self):
      if not self.has_header(): return None
      return SAMHeader(self._header_text)

SAMStreamOptions = namedtuple('SAMStreamOptions',
   ['buffer_size', # reference dictionary
    'reference'])

class SAMStream(object):
  """minimum_intron_size greater than zero will only show sam entries with introns (junctions)
  minimum_overhang greater than zero will require some minimal edge support to consider an intron (junction)

  :param stream: filehandle to go through
  :param options.buffer_size: buffer_size - read this many bytes at a time
  :type stream: stream
  :type options.buffer_size: 
  """
  def __init__(self,stream,options=None):
    if not options: options = SAMStream.Options()
    self._options = options
    self._gen = SAMGenerator(stream,options=self._options)
    self._stream = iter(self._gen)

  @staticmethod
  def Options(**kwargs):
     """ A method for declaring options for the class"""
     construct = SAMStreamOptions #IMPORTANT!  Set this
     names = construct._fields
     d = {}
     for name in names: d[name] = None #default values
     d['buffer_size'] = 10000000
     for k,v in kwargs.iteritems():
       if k in names: d[k] = v
       else: raise ValueError('Error '+k+' is not a property of these options')
     """Create a set of options based on the inputs"""
     return construct(**d)

  def has_header(self): return self._gen.has_header()

  @property
  def header(self):
      if not self.has_header(): return None
      return SAMHeader(self._gen._header_text)

  @property
  def header_text(self):
     return self._gen._header_text

  def __iter__(self):
    return self

  def next(self):
    r = self.read_entry()
    if not r: raise StopIteration
    else: return r

  def read_entry(self):
    try: r = self._stream.next()
    except StopIteration: r = None
    return r

def is_junction_line(line,minlen=68,minoverhang=0):
  prog = re.compile('([0-9]+)([NMX=])')
  f = line.rstrip().split("\t")
  v = prog.findall(f[5])
  #get the indecies of introns
  ns = [i for i in range(0,len(v)) if v[i][1]=='N' and int(v[i][0]) >= minlen]
  if len(ns) == 0: return False
  if minoverhang==0: return True
  good_enough = False
  for intron_index in ns:
    left = sum([int(x[0]) for x in v[0:intron_index] if x[1] != 'N'])
    right = sum([int(x[0]) for x in v[intron_index+1:] if x[1] != 'N'])
    worst = min(left,right)
    if worst >= minoverhang: good_enough = True
  if good_enough: return True
  return False


def is_header(line):
  """true if we are in a header"""
  if re.match('^@',line):
    f = line.rstrip().split("\t")
    if(len(f) > 9):
      return False
    return True
  return False

def sort_header(header_text):
  """sort the chromosomes in a header text"""
  lines = header_text.rstrip().split("\n")
  rlens = {}
  for ln in lines:
    m = re.match('@SQ\tSN:(\S+)\tLN:(\S+)',ln)
    if m:
      rlens[m.group(1)] = m.group(2)
  output = ''
  done_lens = False
  for ln in lines:
    if re.match('@SQ\tSN:',ln):
      if not done_lens:
        done_lens = True
        for chr in sorted(rlens.keys()):
          output += "@SQ\tSN:"+chr+"\tLN:"+str(rlens[chr])+"\n"
    else:
      output += ln.rstrip("\n")+"\n"
  return output

def check_flag(flag,inbit):
  """Check a flag is true or false"""
  if flag & inbit: return True
  return False

