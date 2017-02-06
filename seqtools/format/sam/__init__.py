"""Classes to work with sam and files"""
import sys, re, os, itertools
from collections import namedtuple

import seqtools.align
from seqtools.sequence import rc
from seqtools.range import GenomicRange
import seqtools.stream
_sam_cigar_target_add = re.compile('[M=XDN]$')


"""SAM line options that are not absolutely necessary for a sam line"""
SAMOptions = namedtuple('SAMOptions',
   ['reference', # reference dictionary
    'reference_lengths', # lengths of chromosomes in reference dictionary
    'sam_header',
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
    'tags'])

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
    #self.entries = self._parse_sam_line()
    self._alignment_ranges = None
    self._entries = None
    if self.is_aligned():
       super(SAM,self).__init__(options)
    return

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
    #rname
    rname = None
    if f[2] != '*': rname = f[2]
    #pos
    pos = int(f[3])
    if pos == 0: pos = None 
    #mapq
    mapq = int(f[4])
    if mapq == 255: mapq = None 
    #CIGAR
    cigar = None
    if f[5] != '*': cigar = f[5]
    #rnext
    rnext = None
    if f[6] != '*': rnext = f[6]
    #pnext
    pnext = int(f[7])
    if pnext == 0: pnext = None 
    #tlen
    tlen = int(f[8])
    if tlen == 0: tlen = None 
    #seq
    seq = None
    if f[9] != '*': seq = f[9]
    #qual
    qual = None
    if f[10] != '*': qual = f[10]
    tags = None
    if len(f) > 11:
       tags = f[11:]

    return  SAMFields(
       f[0], #qname
       int(f[1]), #flag
       rname, #rname
       pos, #pos
       mapq, #mapq
       cigar, #cigar
       rnext, #rnext
       pnext, #pnext
       tlen, #tlen
       seq, #seq
       qual, #qual
       tags)

  def __str__(self):
    """The only way to set a SAM line is with a SAM line so we have this"""
    return self._line

  @property
  def cigar(self):
     """cache this one to speed things up a bit"""
     if self._cigar: return self._cigar
     self._cigar = [[int(m[0]),m[1]] for m in re.findall('([0-9]+)([MIDNSHP=X]+)',self.entries.cigar)]
     return self._cigar

  @property
  def tags(self):
     """Access the auxillary data here"""
     if self._tags: return self._tags
     tags = {}
     if not tags: return {}
     for m in [[y.group(1),y.group(2),y.group(3)] for y in [re.match('([^:]{2,2}):([^:]):(.+)$',x) for x in self.entries.tags]]:
        if m[1] == 'i': m[2] = int(m[2])
        elif m[1] == 'f': m[2] = float(m[2])
        tags[m[0]] = {'type':m[1],'value':m[2]}
     self._tags = tags
     return self._tags

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
    if self._reference_lengths:
      if self.entries.rname in self._reference_lengths:
        return self._reference_lengths[self.entries.rname]
    elif self._reference:
      return len(self._reference[self.entries.rname])
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
    if self.check_flag(0x10): return self.entries.qual[::-1]
    return self.entries.qual

  @property
  def query_sequence_length(self):
    """ does not include hard clipped"""
    if seq: return len(self.entries.seq)
    if not self.entries.cigar:
       raise ValueError('Cannot give a query length if no cigar and no query sequence are present')
    return sum([x[0] for x in self.entries.cigar if re.match('[MIS=X]',x[1])])

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
    return sum([x[0] for x in self.entries.get_cigar() if re.match('[HMIS=X]',x[1])])

  @property
  def actual_original_query_range(self):
    """ This accounts for hard clipped bases 
    and a query sequence that hasnt been reverse complemented

    :return: the range covered on the original query sequence
    :rtype: GenomicRange
    """
    l = self.original_query_sequence_length
    a = self.get_alignment_ranges()
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
    cig = [x[:] for x in self.cigar]
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
    tlen = sum([x[0] for x in self.cigar if _sam_cigar_target_add.match(x[1])])
    self._target_range = GenomicRange(self.value('rname'),self.value('pos'),self.value('pos')+tlen-1)
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

class SAMHeader:
  """class to retain information about a SAMheader and access data from it"""
  def __init__(self,header_text):
    self._text = header_text
    self.tags = []
    self._sequence_lengths = {}
    for line in self._text.split("\n"):
      if len(line) == 0: continue
      tag = line[1:3]
      rem = line[4:]
      self.tags.append({'tag':tag,'info':{}})
      for c in [{'field':x[0:2],'value':x[3:]} for x in rem.split("\t")]:
        self.tags[-1]['info'][c['field']] = c['value']
    for v in [x['info'] for x in self.tags if x['tag'] == 'SQ']:
      self._sequence_lengths[v['SN']] = int(v['LN'])
    return

  def __str__(self): return self._text.rstrip()

  def text(self): return self._text.rstrip()+"\n"

  @property
  def sequence_names(self):
    return self._sequence_lengths.keys()

  def get_sequence_lengths(self):
    """return a dictionary of sequence lengths"""
    return self._sequence_lengths

  def get_sequence_length(self,sname):
    """get a specific sequence length"""
    return self._sequence_lengths[sname]


SAMGeneratorOptions = namedtuple('SAMGeneratorOptions',
   ['buffer_size'])

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
      d['buffer_size'] = 1000000
      for k,v in kwargs.iteritems():
         if k in names: d[k] = v
         else: raise ValueError('Error '+k+' is not a property of these options')
      """Create a set of options based on the inputs"""
      return construct(**d)

   def __iter__(self):
      return itertools.imap(SAM,self._gen())

   def has_header(self):
      if len(self._header_text) > 0: return True
      return False

   @property
   def header(self):
      if not self.has_header(): return None
      return SAMHeader(self._header_text)

"""SAM line options that are not absolutely necessary for a sam line"""
SAMStreamOptions = namedtuple('SAMStreamOptions',
   ['reference' # reference dictionary
               ])

class SAMStream:
  """minimum_intron_size greater than zero will only show sam entries with introns (junctions)
  minimum_overhang greater than zero will require some minimal edge support to consider an intron (junction)

  :param fh: filehandle to go through
  :param options.reference: dictionary of reference sequences
  :type fh: stream
  :type options.reference: dict()
  """
  def __init__(self,stream,options=None):
    if not options: options = SAMStream.Options()
    self.previous_line = None
    self.in_header = True
    self._reference = reference
    self._header = None
    self.header_text = ''
    self._stream = stream
    self.assign_handle(stream)
    self.get_header()

  @staticmethod
  def Options(**kwargs):
     """ A method for declaring options for the class"""
     construct = SAMStreamOptions #IMPORTANT!  Set this
     names = construct._fields
     d = {}
     for name in names: d[name] = None #default values
     for k,v in kwargs.iteritems():
       if k in names: d[k] = v
       else: raise ValueError('Error '+k+' is not a property of these options')
     """Create a set of options based on the inputs"""
     return construct(**d)


  # return a string that is the header
  def get_header(self):
    """Return the object representing the header"""
    if not self._header:
      self._header = SAMHeader(self.header_text)
      return self._header
    return self._header


  def set_junction_only(self,mybool=True):
    self.junction_only = mybool

  def assign_handle(self,fh):
    if self.in_header:
      while True:
        self.previous_line = fh.readline()
        if is_header(self.previous_line):
          self.header_text += self.previous_line
          #self.header.append(self.previous_line)
        else:
          self.in_header = False
          self.previous_line = self.previous_line
          break
      # make sure our first line is
      if self.junction_only:
        while True:
          if not self.previous_line: break
          if is_junction_line(self.previous_line,self.minimum_intron_size,self.minimum_overhang): break
          self.previous_line = self.fh.readline()

  def __iter__(self):
    return self

  def next(self):
    r = self.read_entry()
    if not r:
      raise StopIteration
    else:
      return r

  def read_entry(self):
    if not self.previous_line: return False
    out = self.previous_line
    self.previous_line = self.fh.readline()
    if self.junction_only:
      while True:
        if not self.previous_line: break
        if is_junction_line(self.previous_line,self.minimum_intron_size,self.minimum_overhang): break
        self.previous_line = self.fh.readline()
    if out:
      s = SAM(out,reference=self._reference)
      s.get_range()
      return s
    return None

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

