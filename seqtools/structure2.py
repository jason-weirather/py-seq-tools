"""this collection of classes helps us operate on mappings of transcripts"""
import sys, random, string, uuid, pickle, zlib, base64
from seqtools.range import GenomicRange, ranges_to_coverage, merge_ranges
from seqtools.sequence import rc
import seqtools.graph
from math import sqrt

class Transcript:
  """Class to describe the mapping of a basic transcript"""
  def __init__(self):
    self._exons = []
    self._junctions = []
    self._direction = None
    self._transcript_name = None
    self._gene_name = None
    self._name = None # for a single name
    self._range = None # set if not chimeric
    self._id = str(uuid.uuid4())
    self._payload = []
    self._sequence = None

  # _initialize is a dummy function that is run when methods are accessed
  # This allows us to hold of running the time consuming initialization of
  # the GPD until its actually accessed and we can defer this burden more
  # easily in multiprocessing
  def _initialize(self): return  

  @property
  def exons(self):
    """Maybe the most core property of the transcript are the exon defintions.  This is an array of exons."""
    self._initialize()
    return self._exons

  @property
  def junctions(self):
    """Can be inferred from the exons, this is an array of junctions"""
    self._initialize()
    return self._junctions

  def validate(self):
    """be certain the scructure is a transcriptome

    :return: true if exon order is compatible with a transcriptome
    :rtype: list
    """
    self._initialize()
    # check the structure
    prev = None
    for exon in self.exons:
      rng = exon.rng
      if prev:
        if rng.start > prev.end and rng.end >= rng.start:
          continue
        else:
          return False
      prev = rng
    return True

  def copy(self):
    """A copy of the transcript

    :return: transcript copy
    :rtype: Transcript
    """
    self._initialize()
    tx_str = self.dump_serialized() 
    tx = Transcript()
    tx.load_serialized(tx_str)
    return tx

  def subset(self,start,finish):
    """Make a trimmed transcript
       Pre: Start base index 0
       Post: Finish base index 1

    :param start: 0-index start
    :type start: int
    :param end:
    :type end: int
    :return: subset transcript
    :rtype: Transcript
    """
    self._initialize()
    # construct a new transcript
    #print str(start)+' to '+str(finish)
    tx = self.copy()
    keep_ranges = []
    index = 0
    z = 0
    for exon in tx.exons:
      z+=1
      original_rng = exon.rng
      rng = exon.rng.copy()
      done = False;
      #print 'exon length '+str(rng.length())
      if start >= index and start < index+original_rng.length(): # we are in this one
        rng.start = original_rng.start+(start-index) # fix the start
        #print 'fixstart '+str(original_rng.start)+' to '+str(rng.start)
      if finish > index and finish <= index+original_rng.length():
        rng.end = original_rng.start+(finish-index)-1
        done = True
        #print 'fixend '+str(original_rng.end)+' to '+str(rng.end)
 
      if finish <= index+original_rng.length(): # we are in the last exon we need
        index+= original_rng.length()
        keep_ranges.append(rng)
        break
      if index+original_rng.length() < start: # we don't need any bases from this
        index += original_rng.length()
        continue # we don't use this exon
      keep_ranges.append(rng)
      index += original_rng.length()
      if index > finish: break
      if done: break
    tx.set_exons_and_junctions_from_ranges(keep_ranges)
    #print 'ranges:'
    #for rng in keep_ranges:
    #  print rng
    #  print rng.length()
    return tx

  def dump_serialized(self):
    """Generate a string representation of the transcript

    :return: serialized_object
    :rtype: string
    """
    self._initialize()
    ln = self.get_fake_gpd_line()
    return base64.b64encode(zlib.compress(pickle.dumps([ln,self._direction,self._transcript_name,self._gene_name,\
                         self._range,self._id,self._payload,\
                         self._sequence])))
  def load_serialized(self,instr):
    """Load a serialized object string into the object

    :param instr: The serialized string
    :type instr: string
    """
    self._initialize()
    vals = pickle.loads(zlib.decompress(base64.b64decode(instr)))
    import seqtools.format.GPD as inGPD
    gpd = inGPD.GPD(vals[0])
    self.exons = gpd.exons
    self.junctions = gpd.junctions
    self._direction = vals[1]
    self._transcript_name = vals[2]
    self._gene_name = vals[3]
    self._range = vals[4] # set if not chimeric
    self._id = vals[5]
    self._payload = vals[6]
    self._sequence = vals[7]

  def get_junction_string(self):
    """Make a string representation of all the junctions

    :return: junctions as a string
    :rtype: string
    """
    self._initialize()
    if len(self.exons) < 2: return None
    return ",".join([x.get_string() for x in self.junctions])

  def set_payload(self,val):
    """Set a payload for this object

    :param val: payload to be stored
    :type val: Anything that can be put in a list
    """
    self._initialize()
    self._payload = [val]
  def get_payload(self):
    """Get the payload currently being stored

    :return: payload
    :rtype: anything that can be stored in a list
    """
    self._initialize()
    return self._payload[0]

  def get_id(self):
    """Return a unique id created for this transcript when it was made

    :return: uuid4 id as a string
    :rtype: string
    """
    return self._id

  def avg_mutual_coverage(self,gpd):
    """get the coverage fraction of each transcript then return the geometric mean

    :param gpd: Another transcript
    :type gpd: Transcript
    :return: avg_coverage
    :rtype: float
    """
    ov = self.overlap_size(gpd)
    if ov == 0: return 0
    xfrac = float(ov) / float(self.get_length())
    yfrac = float(ov) / float(gpd.get_length())
    return sqrt(xfrac*yfrac)

  # Post: Return the number of overlapping base pairs
  #       between self and tx2 transcript
  def overlap_size(self,tx2):
    """Return the number of overlapping base pairs between two transcripts

    :param tx2: Another transcript
    :type tx2: Transcript
    :return: overlap size in base pairs
    :rtype: int
    """
    self._initialize()
    total = 0
    for e1 in [x.get_range() for x in self.exons]:
      for e2 in [x.get_range() for x in tx2.exons]:
        total += e1.overlap_size(e2)
    return total

  def get_exon_count(self):
    """Count the exons in the transcript

    :return: exon count
    :rtype: int
    """
    self._initialize()
    return len(self.exons)

  def set_range(self):
    """Use the exons that are already present to set the range.  In a refactor this seems like it should go away or become private."""
    self._initialize()
    if len(self.exons) == 0: return None # its ... nothing
    chrs = list(set([x.rng.chr for x in self.exons]))
    if len(chrs) > 1: return None # its chimeric
    self._range = GenomicRange(chrs[0],self.exons[0].rng.start,self.exons[-1].rng.end)

  def get_range(self):
    """Get the range from the leftmost exon to the rightmost

    :return: total range
    :rtype: GenomicRange
    """
    self._initialize()
    if self._range:
      return self._range
    return GenomicRange(self.exons[0].get_range().chr,self.exons[0].get_range().start,self.exons[-1].get_range().end)
  def union(self,tx2): # keep direction and name of self
    """Find the union, or perhaps intersection is a better word for it, for two transcripts.  This makes a new transcript.

    :param tx2: transcript 2
    :type tx2: Transcript
    :return: overlapping portion of the transcripts
    :rtype: Transcript
    """
    self._initialize()
    all = []
    for rng1 in [x.rng for x in self.exons]:
      for rng2 in [y.rng for y in tx2.exons]:
        u = rng1.union(rng2)
        if u: all.append(u)
    if len(all) == 0: return None
    rngs = merge_ranges(all)
    tx = Transcript()
    tx.set_exons_and_junctions_from_ranges(rngs)
    tx._direction = self._direction
    tx._transcript_name = self._transcript_name
    tx._gene_name = self._gene_name
    return tx

  def smooth_gaps(self,min_intron):
    """any gaps smaller than min_intron are joined, andreturns a new transcript with gaps smoothed

    :param min_intron: the smallest an intron can be, smaller gaps will be sealed
    :type min_intron: int
    :return: a mapping with small gaps closed
    :rtype: Transcript
    """
    self._initialize()
    tx = Transcript()
    rngs = [self.exons[0].rng.copy()]
    for i in range(len(self.exons)-1):
      dist = self.exons[i+1].rng.start - rngs[-1].end-1
      if dist >= min_intron:
        rngs.append(self.exons[i+1].rng.copy())
      else:
        rngs[-1].end = self.exons[i+1].rng.end
    tx.set_exons_and_junctions_from_ranges(rngs)
    tx._direction = self._direction
    tx._transcript_name = self._transcript_name
    tx._gene_name = self._gene_name
    return tx

  def set_exons_and_junctions_from_ranges(self,rngs):
    """set all exons and subsequestly juntions from these exon ranges;
       does not set direction of transcript;
       ranges need to be ordered in target order left to right

    This is a core feature for setting up a transcript.

    :param rngs: A list of ranges ordered left to right
    :type rngs: GenomicRange[]

    """
    self._initialize()
    self.exons = []
    self.junctions = []
    for e in rngs:
      ex = Exon(GenomicRange(e.chr,e.start,e.end))
      self.exons.append(ex)
    self.exons[0].set_is_leftmost()
    self.exons[-1].set_is_rightmost()
    for i in range(0,len(self.exons)-1):
      # make a junction
      jx = Junction(GenomicRange(self.exons[i].rng.chr,\
                                 self.exons[i].rng.end,\
                                 self.exons[i].rng.end),\
                    GenomicRange(self.exons[i+1].rng.chr,\
                                 self.exons[i].rng.start,\
                                 self.exons[i+1].rng.start))
      jx.set_exon_left(self.exons[i])
      jx.set_exon_right(self.exons[i+1])
      self.junctions.append(jx)
      self.set_range()
    return

  def get_length(self):
    """Return the length of the transcript in bp. Its the sum of the exons

    :return: length
    :rtype: int
    """
    self._initialize()
    return sum([x.get_length() for x in self.exons])

  def set_strand(self,dir):
    """Set the strand (direction)

    :param dir: direction + or -
    :type dir: char
    """
    self._initialize()
    self._direction = dir
  def get_strand(self):
    """Get the strand

    :return: direction + or -
    :rtype: char
    """
    self._initialize()
    return self._direction

  def get_chrom(self):
    """the reference chromosome. greedy return the first chromosome in exon array

    :return: chromosome
    :rtype: string
    """
    self._initialize()
    if len(self.exons)==0: 
      sys.stderr.write("WARNING can't return chromsome with nothing here\n")
      return None
    return self.exons[0].get_range().chr

  def get_sequence(self,ref_dict=None):
    """A strcutre is defined so get,
    if the sequence is not already there, get the sequence from the reference

    :param ref_dict: reference dictionary (only necessary if sequence has not been set already)
    :type ref_dict: dict()
    """
    self._initialize()
    if self._sequence: return self._sequence
    if not ref_dict:
      sys.stderr.write("ERROR: sequence is not defined and reference is undefined\n")
      sys.exit()
    self.set_sequence(ref_dict)
    return self._sequence

  def set_sequence(self,ref_dict):
    """use the reference dictionary to set the transcript's sequence

    :param ref_dict: reference dictionary
    :type ref_dict: dict()
    """
    self._initialize()
    strand = '+'
    if not self._direction:
      sys.stderr.write("WARNING: no strand information for the transcript\n")
    if self._direction: strand = self._direction
    chr = self.get_chrom()
    seq = ''
    for e in [x.get_range() for x in self.exons]:
      seq += ref_dict[chr][e.start-1:e.end]
    if strand == '-':  seq = rc(seq)
    self._sequence = seq.upper()

  def get_gpd_line(self,transcript_name=None,gene_name=None,strand=None):
    """Get the genpred format string representation of the mapping

    :param transcript_name:
    :param gene_name:
    :param strand:
    :type transcript_name: string
    :type gene_name: string
    :type strand: string
    :return: GPD line
    :rtype: string
    """
    self._initialize()
    tname = self._transcript_name
    gname = self._gene_name
    dir = self._direction
    # check for if we just have a single name
    if not tname and not gname:
      if self._name:
        tname = self._name
        gname = self._name
    if not tname: tname = transcript_name
    if not gname: gname = gene_name
    if not dir: dir = strand
    if not tname or not gname or strand:
      sys.stderr.write("ERROR:  transcript name and gene name and direction must be set to output a gpd line or use get_fake_gpd_line()\n")
    out = ''
    out += tname + "\t"
    out += gname + "\t"
    out += self.exons[0].rng.chr + "\t"
    out += dir + "\t"
    out += str(self.exons[0].rng.start-1) + "\t"
    out += str(self.exons[-1].rng.end) + "\t"
    out += str(self.exons[0].rng.start-1) + "\t"
    out += str(self.exons[-1].rng.end) + "\t"
    out += str(len(self.exons)) + "\t"
    out += str(','.join([str(x.rng.start-1) for x in self.exons]))+','+"\t"
    out += str(','.join([str(x.rng.end) for x in self.exons]))+','
    return out

  def set_gene_name(self,name):
    """assign a gene name

    :param name: name
    :type name: string
    """
    self._initialize()
    self._gene_name = name
  def get_gene_name(self):
    """retrieve the gene name

    :return: gene name
    :rtype: string
    """
    self._initialize()
    return self._gene_name
  def set_transcript_name(self,name):
    """assign a transcript name

    :param name: name
    :type name: string
    """
    self._initialize()
    self._transcript_name = name
  def get_transcript_name(self):
    """retrieve the transcript name

    :return: transcript name
    :rtype: string
    """
    self._initialize()
    return self._transcript_name

  def get_fake_psl_line(self,ref):
    """Convert a mapping to a fake PSL line

    :param ref: reference genome dictionary
    :type ref: dict()
    :return: psl line
    :rtype: string
    """
    self._initialize()
    e = self
    mylen = 0
    matches = 0
    qstartslist = []
    for exon in self.exons:
      mylen = exon.rng.length()
      matches += mylen
      qstartslist.append(matches-mylen)
    qstarts = ','.join([str(x) for x in qstartslist])+','
    oline =  str(matches)+"\t" # 1
    oline += "0\t" # 2
    oline += "0\t" # 3
    oline += "0\t" # 4
    oline += "0\t" # 5
    oline += "0\t" # 6
    oline += "0\t" # 7
    oline += "0\t" # 8
    oline += e.get_strand()+"\t" # 9
    oline += e.get_transcript_name()+"\t" # 10
    oline += str(matches)+"\t" # 11
    oline += "0\t" # 12
    oline += str(matches)+"\t" # 13
    oline += e.get_chrom()+"\t" # 14
    oline += str(len(ref[e.get_chrom()]))+"\t" # 15
    oline += str(e.exons[0].rng.start-1)+"\t" # 16
    oline += str(e.exons[-1].rng.end)+"\t" # 17
    oline += str(len(e.exons))+"\t" # 18
    oline += ','.join([str(e.exons[x].rng.end-(e.exons[x].rng.start-1)) for x in range(0,len(e.exons))])+','+"\t" # 19
    oline += qstarts + "\t" # 20
    oline += ','.join([str(x.rng.start-1) for x in e.exons])+',' # 21
    return oline

  def get_fake_gpd_line(self):
    """Convert a mapping to a fake GPD line. not sure why its called fake

    :return: gpd line
    :rtype: string
    """

    self._initialize()
    rlen = 8
    #name = ''.join(random.choice(string.letters+string.digits) for i in range(0,rlen))
    name = str(self.get_id())
    out = ''
    out += name + "\t"
    out += name + "\t"
    out += self.exons[0].rng.chr + "\t"
    out += '+' + "\t"
    out += str(self.exons[0].rng.start-1) + "\t"
    out += str(self.exons[-1].rng.end) + "\t"
    out += str(self.exons[0].rng.start-1) + "\t"
    out += str(self.exons[-1].rng.end) + "\t"
    out += str(len(self.exons)) + "\t"
    out += str(','.join([str(x.rng.start-1) for x in self.exons]))+','+"\t"
    out += str(','.join([str(x.rng.end) for x in self.exons]))+','
    return out

  def get_junctions_string(self):
    """Get a string representation of the junctions.  This is almost identical to a previous function.

    That function is get_junction_string.  A refactor should clear this redundancy.

    :return: string representation of junction
    :rtype: string
    """
    self._initialize()
    return ';'.join([x.get_range_string() for x in self.junctions])

  def junction_overlap(self,tx,tolerance=0):
    """Calculate the junction overlap between two transcripts

    :param tx: Other transcript
    :type tx: Transcript
    :param tolerance: how close to consider two junctions as overlapped (default=0)
    :type tolerance: int
    :return: Junction Overlap Report
    :rtype: Transcript.JunctionOverlap
    """
    self._initialize()
    return Transcript.JunctionOverlap(self,tx,tolerance)

  def exon_overlap(self,tx,multi_minover=10,multi_endfrac=0,multi_midfrac=0.8,single_minover=50,single_frac=0.5,multi_consec=True):
    """Get a report on how mucht the exons overlap

    :param tx:
    :param multi_minover: multi-exons need to overlap by at lest this much to be considered overlapped (default 10)
    :param multi_endfrac: multi-exons need an end fraction coverage of at least this by default (default 0)
    :param multi_midfrac: multi-exons need (default 0.8) mutual coverage for internal exons
    :parma single_minover: single-exons need at least this much shared overlap (default 50)
    :param single_frac: at least this fraction of single exons must overlap (default 0.5)
    :parma multi_consec: exons need to have multiexon consecutive mapping to consider it a match (default True)
    :type tx:
    :type multi_minover: int
    :type multi_endfrac: float
    :type multi_midfrac: float
    :type single_minover: int
    :type single_frac: float
    :type multi_consec: bool
    :return: ExonOverlap report
    :rtype: Transcript.ExonOverlap
    """
    self._initialize()
    return Transcript.ExonOverlap(self,tx,multi_minover,multi_endfrac,multi_midfrac,single_minover,single_frac,multi_consec=multi_consec)

  class ExonOverlap:
    """class to describe exon overlap

    :param tx:
    :param multi_minover: multi-exons need to overlap by at lest this much to be considered overlapped (default 10)
    :param multi_endfrac: multi-exons need an end fraction coverage of at least this by default (default 0)
    :param multi_midfrac: multi-exons need (default 0.8) mutual coverage for internal exons
    :parma single_minover: single-exons need at least this much shared overlap (default 50)
    :param single_frac: at least this fraction of single exons must overlap (default 0.5)
    :parma multi_consec: exons need to have multiexon consecutive mapping to consider it a match (default True)
    :type tx:
    :type multi_minover: int
    :type multi_endfrac: float
    :type multi_midfrac: float
    :type single_minover: int
    :type single_frac: float
    :type multi_consec: bool
    :return: ExonOverlap report
    :rtype: Transcript.ExonOverlap
    """
    def __init__(self1,tx_obj1,tx_obj2,multi_minover=10,multi_endfrac=0,multi_midfrac=0.8,single_minover=50,single_frac=0.5,multi_consec=True):
      self1.tx_obj1 = tx_obj1
      self1.tx_obj2 = tx_obj2
      self1.multi_minover = multi_minover # multi-exon minimum overlap of each exon
      self1.multi_endfrac = multi_endfrac # multi-exon minimum fractional overlap of first or last exon
      self1.multi_midfrac = multi_midfrac # multi-exon minimum fractional overlap of internal exons
      self1.multi_consec = multi_consec # require consecutive exons for exon overlap of multi_exon
      self1.single_minover = single_minover # single-exon minimum overlap
      self1.single_frac = single_frac #single-exon minimum overlap
      self1.overs = [] # set by calculate_overlap()
      #self1.dif1 = []
      #self1.dif2 = []
      self1.calculate_overlap()
      if len(self1.overs) == 0: return None# nothing to analyze
      if self1.tx_obj1.get_exon_count() > 1 and self1.tx_obj1.get_exon_count() > 1 \
         and self1.multi_consec and len(self1.overs) < 2: 
        return None #not enough to consider multi exon transcript overlap
      self1.analyze_overs()
      if self1.tx_obj1.get_exon_count() > 1 and self1.tx_obj1.get_exon_count() > 1 \
        and self1.multi_consec and (min(self1.dif1) != 1 or min(self1.dif2) !=1): 
        return None #not enough to consider multi exon transcript overlap
      
    def __nonzero__(self1):
      if len(self1.overs) > 0: return True
      return False
  
    def match_exon_count(self1):
      """Total number of exons that overlap

      :return: matched exon count
      :rtype: int
      """
      return len(self1.overs)

    def consecutive_exon_count(self1):
      """Best number of consecutive exons that overlap

      :return: matched consecutive exon count
      :rtype: int
      """
      best = 1
      consec = 1
      for i in range(0,len(self1.dif1)):
        if self1.dif1[i] == 1 and self1.dif2[i] == 1:
          consec += 1
        else: 
          consec = 1
        if consec > best:
          best = consec
      return best

    def is_subset(self1):
      """ Return value if tx_obj2 is a complete subset of tx_obj1 or 
          tx_obj1 is a complete subset of tx_obj2

          Values are:

          * Return 1: Full overlap (mutual subests) 
          * Return 2: two is a subset of one
          * Return 3: one is a subset of two
          * Return False if neither is a subset of the other

      :return: subset value
      :rtype: int
      """
      if len(self1.overs) == 0: return False
      if len(self1.dif1) > 0: # make sure they are consecutive if more than one
        if max(self1.dif1) != 1 or max(self1.dif2) != 1: return False
      onecov = self1.start1 and self1.end1
      twocov = self1.start2 and self1.end2
      if onecov and twocov:
        return 1
      elif twocov: return 2
      elif onecov: return 3
      return False

    def is_full_overlap(self1):
      """true if they are a full overlap

      :return: is full overlap
      :rtype: bool
      """
      if len(self1.overs) == 0: return False
      if len(self1.dif1) > 0:
        if max(self1.dif1) != 1 or max(self1.dif2) != 1: return False
      if self1.start1 and self1.end1 and self1.start2 and self1.end2:
        return True
      return False

    def is_compatible(self1):
      """ Return True if the transcripts can be combined together 

      :return: can be combined together
      :rtype: bool
      """
      if len(self1.overs) == 0: return False
      if len(self1.dif1) > 0:
        if max(self1.dif1) != 1 or max(self1.dif2) != 1: return False
      # If we are still here it is a single run
      if (self1.start1 or self1.start2) and (self1.end1 or self1.end2):
        return True
      return False

    def analyze_overs(self1):
      """A helper function that prepares overlap and consecutive matches data"""
      #check for full overlap first
      self1.dif1 = [self1.overs[i][0]-self1.overs[i-1][0] for i in range(1,len(self1.overs))]
      self1.dif2 = [self1.overs[i][1]-self1.overs[i-1][1] for i in range(1,len(self1.overs))]
      #see if it starts and ends on first or last junction
      self1.start1 = self1.overs[0][0] == 0
      self1.start2 = self1.overs[0][1] == 0
      self1.end1 = self1.overs[-1][0] == len(self1.tx_obj1.exons)-1
      self1.end2 = self1.overs[-1][1] == len(self1.tx_obj2.exons)-1
      return

    def calculate_overlap(self1):
      """Create the array that describes how junctions overlap"""
      overs = []
      if not self1.tx_obj1.get_range().overlaps(self1.tx_obj2.get_range()): return # if they dont overlap wont find anything
      for i in range(0,len(self1.tx_obj1.exons)):
        for j in range(0,len(self1.tx_obj2.exons)):
          osize = self1.tx_obj1.exons[i].rng.overlap_size(self1.tx_obj2.exons[j].rng)
          ofrac = 0
          if osize > 0:
            ofrac = min(float(osize)/float(self1.tx_obj1.exons[i].rng.length())\
                       ,float(osize)/float(self1.tx_obj2.exons[j].rng.length()))

          if self1.tx_obj1.get_exon_count() == 1 or self1.tx_obj2.get_exon_count == 1:
            # use single exon rules
            if osize >= self1.single_minover and ofrac >= self1.single_frac:
              #print 'single exon match'
              overs.append([i,j])
          else: # for multi exons
            if i == 0 or j == 0 or i == len(self1.tx_obj1.exons)-1 or j == len(self1.tx_obj2.exons)-1:
              #its on an end
              if osize >= self1.multi_minover and ofrac >= self1.multi_endfrac:
                #print 'end exon match'
                overs.append([i,j])
            #else its a middle
            elif osize >= self1.multi_minover and ofrac >= self1.multi_midfrac:
              #print 'mid exon match'
              overs.append([i,j])
      #print overs
      self1.overs = overs

  class JunctionOverlap:
    """Class for describing the overlap of junctions between transcripts
 
    This should probably be not a child.

    :param tx_obj1: transcript1
    :param tx_obj2: transcript2
    :param tolerance: how far before its no longer a matched junction
    :type tx_obj1: Transcript
    :type tx_obj2: Transcript
    :type tolerance: int
    """
    def __init__(self1,tx_obj1,tx_obj2,tolerance=0):
      self1.tx_obj1 = tx_obj1
      self1.tx_obj2 = tx_obj2
      self1.tolerance = tolerance
      self1.overs = [] # gets set by calculate_overlap()
      self1.calculate_overlap()
      if len(self1.overs) == 0: return None# nothing to analyze
      self1.analyze_overs()

    def __nonzero__(self1):
      if len(self1.overs) > 0: return True
      return False

    def match_junction_count(self1):
      return len(self1.overs)

    def is_subset(self1):
      """Return value if tx_obj2 is a complete subset of tx_obj1 or tx_obj1 is a complete subset of tx_obj2

      values:

      * Return 1: Full overlap (mutual subests) 
      * Return 2: two is a subset of one
      * Return 3: one is a subset of two
      * Return False if neither is a subset of the other
      """
      if len(self1.overs) == 0: return False
      if len(self1.dif1) > 0: # make sure they are consecutive if more than one
        if max(self1.dif1) != 1 or max(self1.dif2) != 1: return False
      onecov = self1.start1 and self1.end1
      twocov = self1.start2 and self1.end2
      if onecov and twocov:
        return 1
      elif twocov: return 2
      elif onecov: return 3
      return False

    def is_full_overlap(self1):
      """True if its a full overlap

      :return: True if its a full overlap
      :rtype: bool
      """
      if len(self1.overs) == 0: return False
      if len(self1.dif1) > 0:
        if max(self1.dif1) != 1 or max(self1.dif2) != 1: return False
      if self1.start1 and self1.end1 and self1.start2 and self1.end2:
        return True
      return False

    def is_compatible(self1):
      """Return True if the transcripts can be combined together

      :return: True if we can combine
      :rtype: bool
      """
      if len(self1.overs) == 0: return False
      if len(self1.dif1) > 0:
        if max(self1.dif1) != 1 or max(self1.dif2) != 1: return False
      # If we are still here it is a single run
      if (self1.start1 or self1.start2) and (self1.end1 or self1.end2):
        return True
      return False

    def analyze_overs(self1):
      """A helper function to prepare values describing overlaps"""
      #check for full overlap first
      self1.dif1 = [self1.overs[i][0]-self1.overs[i-1][0] for i in range(1,len(self1.overs))]
      self1.dif2 = [self1.overs[i][1]-self1.overs[i-1][1] for i in range(1,len(self1.overs))]
      #see if it starts and ends on first or last junction
      self1.start1 = self1.overs[0][0] == 0
      self1.start2 = self1.overs[0][1] == 0
      self1.end1 = self1.overs[-1][0] == len(self1.tx_obj1.junctions)-1
      self1.end2 = self1.overs[-1][1] == len(self1.tx_obj2.junctions)-1
      return

    def calculate_overlap(self1):
      """Create the array that describes how junctions overlap"""
      overs = []
      if not self1.tx_obj1.get_range().overlaps(self1.tx_obj2.get_range()): return # if they dont overlap wont find anything
      for i in range(0,len(self1.tx_obj1.junctions)):
        for j in range(0,len(self1.tx_obj2.junctions)):
          if self1.tx_obj1.junctions[i].overlaps(self1.tx_obj2.junctions[j],self1.tolerance):
            overs.append([i,j])
      self1.overs = overs

class Junction:
  """ class to describe a junction

  :param rng_left: left side of junction
  :param rng_right: right side of junction
  :type rng_left: GenomicRange
  :type rng_right: GenomicRange
  """
  def __init__(self,rng_left=None,rng_right=None):
    self.left = rng_left
    self.right = rng_right
    self.left_exon = None
    self.right_exon = None
  def dump_serialized(self):
    """Get string representation of the junction

    :return: serialized object
    :rtype: string
    """
    return pickle.dumps(self)
  def load_serialized(self,instr):
    """load the string

    :param instr:
    :type instr: string
    """
    self = pickle.loads(instr)
  def get_string(self):
    """A string representation of the junction

    :return: string represnetation
    :rtype: string
    """
    return self.left.chr+':'+str(self.left.end)+'-'+self.right.chr+':'+str(self.right.start)
  def get_left_exon(self):
    """ get the exon to the left of the junction

    :return: left exon
    :rtype: Exon
    """
    return self.left_exon
  def get_right_exon(self):
    """ get the exon to the right of the junction

    :return: right exon
    :rtype: Exon or GenomicRange
    """
    return self.right_exon
  def get_range_string(self):
    """Another string representation of the junction.  these may be redundant."""
    return self.left.chr+":"+str(self.left.end)+'/'+self.right.chr+":"+str(self.right.start)
  def set_left(self,rng):
    """ Assign the leftmost range"""
    self.left = rng
  def set_right(self,rng):
    """ Assign the right most range"""
    self.right = rng
  def equals(self,junc):
    """test equality with another junction"""
    if self.left.equals(junc.left): return False
    if self.right.equals(junc.right): return False
    return True
  def overlaps(self,junc,tolerance=0):
    """see if junction overlaps with tolerance""" 
    if not self.left.overlaps_with_padding(junc.left,tolerance): return False
    if not self.right.overlaps_with_padding(junc.right,tolerance): return False
    return True
  def cmp(self,junc,tolerance=0):
    """ output comparison and allow for tolerance if desired

    * -1 if junc comes before self
    * 1 if junc comes after self
    * 0 if overlaps
    * 2 if else

    :param junc:
    :param tolerance: optional search space (default=0, no tolerance)
    :type junc: Junction
    :type tolerance: int
    :return: value of comparison
    :rtype: int
    """
    if self.overlaps(junc,tolerance):
      return 0 #equal
    if self.left.chr == junc.right.chr:
      if self.left.start > junc.right.start:
        return -1 #comes before
    if self.right.chr == junc.left.chr:
      if self.right.start < junc.right.start:
        return 1 #comes after
    return 2

  def set_exon_left(self,ex):
    """assign the left exon"""
    self.left_exon = ex
    ex.right_junc = self
  def set_exon_right(self,ex):
    """assign the right exon"""
    self.right_exon = ex
    ex.left_junc = self

class Exon:
  """class to describe an exon

  :param rng:
  :type rng: GenomicRange
  """
  def __init__(self,rng=None):
    self.rng = rng
    self.left_junc = None
    self.right_junc = None
    self._is_leftmost = False #bool is it a start or end
    self._is_rightmost = False
  def dump_serialized(self):
    return pickle.dumps(self)
  def load_serialized(self,instr):
    self = pickle.loads(instr)
  def get_range(self):
    return self.rng
  def get_length(self):
    return self.rng.length()
  def set_left_junc(self,junc):
    self.left_junc = junc
    junc.set_right_exon = self
  def set_right_junc(self,junc):
    self.right_junc = junc
    junc.set_left_exon = self
  def set_is_leftmost(self,boo=True): 
    self._is_leftmost = boo # is leftmost
  def set_is_rightmost(self,boo=True): 
    self._is_rightmost = boo   #is rightmost

class TranscriptLoci:
  """combine together compatible multiple transcript groups to form
   a simpler set of transcripts """
  def __init__(self):
    #self.transcripts = []
    self.merge_rules = TranscriptLociMergeRules('is_any_overlap')
    self.merge_rules.set_juntol(10)
    self.g = seqtools.Graph.Graph()   

  def __str__(self):
    return str(len(self.g.get_nodes()))+ " nodes"  

  def remove_transcript(self,tx_id):
    """Remove a transcript from the locus by its id

    :param tx_id:
    :type tx_id: string
    """
    txs = self.get_transcripts()
    if tx_id not in [x.get_id() for x in txs]:
      return
    tx = [x for x in txs if x.get_id()==tx_id][0]
    for n in [x for x in self.g.get_nodes()]:
      if tx_id not in [y.get_id() for y in n.get_payload()]:
        continue
      n.get_payload().remove(tx)
      if len(n.get_payload())==0:
        self.g.remove_node(n)      
      #sys.stderr.write("\n"+str(n.get_payload())+"\n")
      # 
      #sys.stderr.write("\n"+str(len(n.get_payload()))+"\n")
  def set_merge_rules(self,mr):  
    """Define rules for how to merge transcripts

    :param mr:
    :type mr: TranscriptLociMergeRules
    """
    self.merge_rules = mr

  def get_depth_per_transcript(self,mindepth=1):
    """ using all the transcripts find the depth """
    bedarray = []
    for tx in self.get_transcripts():
      for ex in [x.rng for x in tx.exons]: bedarray.append(ex)
    cov = ranges_to_coverage(bedarray)
    results = {}
    for tx in self.get_transcripts():
      tlen = tx.get_length()
      bcov = []
      for ex in [x.rng for x in tx.exons]:     
        excov = [[x.overlap_size(ex),x.get_payload()] for x in cov]
        for coved in [x for x in excov if x[0] > 0]:
          bcov.append(coved)
      total_base_coverage = sum([x[0]*x[1] for x in bcov])
      average_coverage = float(total_base_coverage)/float(tlen)
      minimum_bases_covered = sum([x[0] for x in bcov if x[1] >= mindepth])
      fraction_covered_at_minimum = float(minimum_bases_covered)/float(tlen)
      res = {'tx':tx,'average_coverage':average_coverage,'fraction_covered':fraction_covered_at_minimum,'mindepth':mindepth,'length_covered':minimum_bases_covered}
      results[tx.get_id()] = res
      #print average_coverage
      #print fraction_covered_at_minimum
      #print tlen
      #tcov = float(bcov)/float(tlen)
      #print tcov
    #for c in cov:
    #  print c
    return results

  def get_range(self):
    """Return the range the transcript loci covers

    :return: range
    :rtype: GenomicRange
    """
    chrs = set([x.get_range().chr for x in self.get_transcripts()])
    if len(chrs) != 1: return None
    start = min([x.get_range().start for x in self.get_transcripts()])
    end = max([x.get_range().end for x in self.get_transcripts()])
    return GenomicRange(list(chrs)[0],start,end)

  def get_transcripts(self):
    """ a list of the transcripts in the locus"""
    txs = []
    for pays in [x.get_payload() for x in self.g.get_nodes()]:
      for pay in pays:
        txs.append(pay)
    return txs

  def partition_loci(self,verbose=False):
    """ break the locus up into unconnected loci

    :return: list of loci
    :rtype: TranscriptLoci[]
    """
    #names = []
    #for entries in [x.get_payload() for x in self.g.get_nodes()]:
    #  for entry in entries:
    #    names.append(entry.get_gene_name())

    #sys.stderr.write('-------partition_loci-----'+"\n")
    #sys.stderr.write(self.g.get_report()+"\n")
    self.g.merge_cycles()
    #sys.stderr.write(self.g.get_report()+"\n")
    gs = self.g.partition_graph(verbose=verbose)
    tls = [] # makea list of transcript loci
    for g in gs:
      tl = TranscriptLoci()
      tl.merge_rules = self.merge_rules
      ns = g.get_nodes()
      for n in [x.get_payload() for x in ns]:
        for tx in n:
          tl.add_transcript(tx)
      if len(tl.g.get_nodes()) > 0:
        tls.append(tl)
    #print '-----------------------' 
    #names = []
    #for tl in tls:
    #  for tx in tl.get_transcripts():
    #    names.append(tx.get_gene_name())
    #for name in sorted(names):
    #  print name
    #print '--------------------------'
    return tls

  def add_transcript(self,tx):
    """Add a transcript to the locus

    :param tx: transcript to add
    :type tx: Transcript
    """
    for y in [x.get_payload() for x in self.g.get_nodes()]:
      if tx.get_id in [z.get_id() for z in y]:
        #if tx.get_id() in [[y.get_id() for y in x.get_payload()] for x in self.g.get_nodes()]:
        sys.stderr.write("WARNING tx is already in graph\n")
        return True
    # transcript isn't part of graph yet
    n = seqtools.Graph.Node([tx])

    other_nodes = self.g.get_nodes()
    self.g.add_node(n)
    # now we need to see if its connected anywhere
    for n2 in other_nodes:
     tx2s = n2.get_payload()
     for tx2 in tx2s:
      # do exon overlap
      er = self.merge_rules.get_exon_rules()
      # if we are doing things by exon
      if (self.merge_rules.get_use_single_exons() and (tx.get_exon_count() == 1 or tx2.get_exon_count() == 1)) or \
         (self.merge_rules.get_use_multi_exons() and (tx.get_exon_count() > 1 and tx2.get_exon_count() > 1)):
        eo = tx.exon_overlap(tx2,multi_minover=er['multi_minover'],multi_endfrac=er['multi_endfrac'],multi_midfrac=er['multi_midfrac'],single_minover=er['single_minover'],single_frac=er['single_frac'])
        if self.merge_rules.get_merge_type() == 'is_compatible':
          if eo.is_compatible():
            self.g.add_edge(seqtools.Graph.Edge(n,n2),verbose=False)
            self.g.add_edge(seqtools.Graph.Edge(n2,n),verbose=False)
        elif self.merge_rules.get_merge_type() == 'is_subset':
          r = eo.is_subset()
          if r == 2 or r == 1:
            self.g.add_edge(seqtools.Graph.Edge(n,n2),verbose=False)
          if r == 3 or r == 1:
            self.g.add_edge(seqtools.Graph.Edge(n2,n),verbose=False)
        elif self.merge_rules.get_merge_type() == 'is_full_overlap':
          if eo.is_full_overlap():
            self.g.add_edge(seqtools.Graph.Edge(n,n2),verbose=False)
            self.g.add_edge(seqtools.Graph.Edge(n2,n),verbose=False)
        elif self.merge_rules.get_merge_type() == 'is_any_overlap':
          if eo.match_exon_count() > 0:
            self.g.add_edge(seqtools.Graph.Edge(n,n2),verbose=False)
            self.g.add_edge(seqtools.Graph.Edge(n2,n),verbose=False)        
            
      if self.merge_rules.get_use_junctions():
        # do junction overlap
        jo = tx.junction_overlap(tx2,self.merge_rules.get_juntol())
        #print jo.match_junction_count()
        if self.merge_rules.get_merge_type() == 'is_compatible':
          if jo.is_compatible():
            self.g.add_edge(seqtools.Graph.Edge(n,n2),verbose=False)
            self.g.add_edge(seqtools.Graph.Edge(n2,n),verbose=False)
        elif self.merge_rules.get_merge_type() == 'is_subset':
          r = jo.is_subset()
          if r == 2 or r == 1:
            self.g.add_edge(seqtools.Graph.Edge(n,n2),verbose=False)
          if r == 3 or r == 1:
            self.g.add_edge(Seqtools.Graph.Edge(n2,n),verbose=False)
        elif self.merge_rules.get_merge_type() == 'is_full_overlap':
          if jo.is_full_overlap():
            self.g.add_edge(seqtools.Graph.Edge(n,n2),verbose=False)
            self.g.add_edge(seqtools.Graph.Edge(n2,n),verbose=False)
        elif self.merge_rules.get_merge_type() == 'is_any_overlap':
          if jo.match_junction_count() > 0:
            self.g.add_edge(seqtools.Graph.Edge(n,n2),verbose=False)
            self.g.add_edge(seqtools.Graph.Edge(n2,n),verbose=False)        
    return True
  #def add_transcript_group(self,txg):
  #  self.transcript_groups.append(txg)      

  #def merge_down_loci(self):
  #  # look at the transcript groups that are currently there
  #  # check for full match
  #  
  #  return

# TranscriptLocus Merge Rules
class TranscriptLociMergeRules:
    """Establish rules up on which to merge loci"""
    def __init__(self,merge_type):
      #Multi-exon rules
      self._junction_tolerance = 10
      self._possible_types = set(['is_subset','is_compatible','is_full_overlap','is_any_overlap'])
      if merge_type not in self._possible_types:
        sys.stderr.write("ERROR: "+merge_type+" is not a known merge type\n")
        sys.exit()
      self._merge_type = merge_type
      self._use_junctions = True
      # exon rules
      self._use_multi_exons = True
      self._use_single_exons = True
      self._multi_minover=10
      self._multi_endfrac=0
      self._multi_midfrac=0.8
      self._single_minover=100
      self._single_frac=0.5
      return
    def get_use_single_exons(self): return self._use_single_exons
    def get_use_multi_exons(self): return self._use_multi_exons
    def get_use_junctions(self): return self._use_junctions
    def get_exon_rules(self):
      return {'multi_minover':self._multi_minover,'multi_endfrac':self._multi_endfrac,'multi_midfrac':self._multi_midfrac,'single_minover':self._single_minover,'single_frac':self._single_frac}
    def get_merge_type(self):
      return self._merge_type
    def set_juntol(self,juntol):
      self._junction_tolerance = juntol
    def get_juntol(self):
      return self._junction_tolerance
    def set_use_junctions(self,boo=True):
      self._use_junctions = boo

class TranscriptGroup:
  """A transcript group is like the fuzzy gpd class we had before"""
  def __init__(self):
    self.junction_groups = [] # These will be more fuzzy defitions
    #self.exons = [] # These will be based on the junctions and individual starts
    self.transcripts = [] # These are the individual transcripts that make up this group
    self._transcript_ids = set()

  def get_transcript(self,exon_bounds='max'):
    """Return a representative transcript object"""
    out = Transcript()
    out.junctions = [x.get_junction() for x in self.junction_groups]
    # check for single exon transcript
    if len(out.junctions) == 0:
      leftcoord = min([x.exons[0].rng.start for x in self.transcripts])
      rightcoord = max([x.exons[-1].rng.end for x in self.transcripts])
      e = Exon(GenomicRange(x.exons[0].rng.chr,leftcoord,rightcoord))
      e.set_is_leftmost()
      e.set_is_rightmost()
      out.exons.append(e)
      return out
    # get internal exons
    self.exons = []
    for i in range(0,len(self.junction_groups)-1):
      j1 = self.junction_groups[i].get_junction()
      j2 = self.junction_groups[i+1].get_junction()
      e = Exon(GenomicRange(j1.right.chr,j1.right.end,j2.left.start))
      e.set_left_junc(j1)
      e.set_right_junc(j2)
      #print str(i)+" to "+str(i+1)
      out.exons.append(e)
    # get left exon
    left_exons = [y for y in [self.transcripts[e[0]].junctions[e[1]].get_left_exon() for e in self.junction_groups[0].evidence] if y]
    if len(left_exons) == 0:
      sys.stderr.write("ERROR no left exon\n")
      sys.exit()
    e_left = Exon(GenomicRange(out.junctions[0].left.chr,\
                               min([x.get_range().start for x in left_exons]),
                               out.junctions[0].left.start))
    e_left.set_right_junc(out.junctions[0])
    out.exons.insert(0,e_left)
    # get right exon
    right_exons = [y for y in [self.transcripts[e[0]].junctions[e[1]].get_right_exon() for e in self.junction_groups[-1].evidence] if y]
    if len(right_exons) == 0:
      sys.stderr.write("ERROR no right exon\n")
      sys.exit()
    e_right = Exon(GenomicRange(out.junctions[-1].right.chr,\
                               out.junctions[-1].right.end,\
                               max([x.get_range().end for x in right_exons])))
    e_right.set_left_junc(out.junctions[-1])
    out.exons.append(e_right)
    return out

  def add_transcript(self,tx,juntol=0,verbose=True):
    if tx.get_id() in self._transcript_ids: return True
    # check existing transcripts for compatability
    for t in self.transcripts:
      ov = t.junction_overlap(tx,juntol)
      if ov:
        if not ov.is_compatible():
          if verbose: sys.stderr.write("transcript is not compatible\n") 
          return False

      else: 
        if verbose: sys.stderr.write("transcript is not overlapped\n")
        return False # if its not overlapped we also can't add
    self.transcripts.append(tx)
    curr_tx = len(self.transcripts)-1
    #print curr_tx
    # see if there is no junctions yet
    if len(self.junction_groups) == 0:
      for i in range(0,len(tx.junctions)):
        jg = TranscriptGroup.JunctionGroup(self)
        jg.add_junction(curr_tx,i,tolerance=juntol)
        self.junction_groups.append(jg)
    else: # there is already a transcript(s) here to work around
      before = []
      middle = []
      after = []
      for j in range(0,len(tx.junctions)):
        # see if its before the existing set
        cmp = self.junction_groups[0].get_junction().cmp(tx.junctions[j])
        if cmp == -1: before.append(j)
        # see if it goes in the existing set
        for k in range(0,len(self.junction_groups)):
          ov = self.junction_groups[k].get_junction().overlaps(tx.junctions[j],tolerance=juntol) #may need to add a tolerance
          if ov: middle.append([j,k])
        # see if it goes after this set
        cmp = self.junction_groups[-1].get_junction().cmp(tx.junctions[j])
        if cmp == 1: after.append(j)
      # add to the middle values before we disrupt indexing
      #print '---'
      #print len(before)
      #print len(middle)
      #print len(after)
      #print '---'
      for v in middle:
        self.junction_groups[v[1]].add_junction(curr_tx,v[0],tolerance=juntol) 
      #add to the beginning and then the end
      for i in reversed(before):
         jg = TranscriptGroup.JunctionGroup(self)
         jg.add_junction(curr_tx,i,tolerance=juntol)
         self.junction_groups.insert(0,jg)        
      for i in after:
         jg = TranscriptGroup.JunctionGroup(self)
         jg.add_junction(curr_tx,i,tolerance=juntol)
         self.junction_groups.append(jg)        
      #if len(tx.junctions)==0:
      #  jg = TranscriptGroup.JunctionGroup(self)
      #  jg.add_junction(curr_tx,i)
      #  self.junctions.append(jg)
    self._transcript_ids.add(tx.get_id())
    return True

  class JunctionGroup:
    """Describe a junction as a group of junctions with options for junction tolerance"""
    def __init__(self1,outer):
      self1.outer = outer
      self1.evidence = [] # array of evidence that is the 
                         # outer.transcript index
                         # outer.trascript.junction index
      self1.representative_junction = None #calculated as needed
    def get_junction(self1): 
      """return the consensus junction"""
      if self1.representative_junction:
        return self1.representative_junction
      left_rngs = []
      right_rngs = []
      for j in [self1.outer.transcripts[x[0]].junctions[x[1]] for x in self1.evidence]:
        left_rngs.append(j.left)
        right_rngs.append(j.right)
      left = _mode([x.end for x in left_rngs])
      right = _mode([x.start for x in right_rngs])
      outj = Junction(GenomicRange(left_rngs[0].chr,left,left),GenomicRange(right_rngs[0].chr,right,right))
      self1.representative_junction = outj
      return outj
    def add_junction(self1,tx_index,junc_index,tolerance=0):
      """add a junction"""
      self1.representative_junction = None
      if len(self1.evidence)==0:
        # go ahead and add it
        #j = self1.outer.transcripts[tx_index].junctions[junc_index]
        self1.evidence.append([tx_index,junc_index])
      else:
        # check it and add it
        if not self1.get_junction().overlaps(self1.outer.transcripts[tx_index].junctions[junc_index],tolerance=tolerance):
          sys.stderr.write("WARNING Unable to add junction JunctionGroup\n"+self1.get_junction().get_range_string()+"\n"+self1.outer.transcripts[tx_index].junctions[junc_index].get_range_string()+"\n")
          return False
        self1.evidence.append([tx_index,junc_index])

def _mode(mylist):
  counts = [mylist.count(x) for x in mylist]
  maxcount = max(counts)
  avg = sum([float(x) for x in mylist])/len(mylist)
  #print counts 
  dist = [abs(float(x)-avg) for x in mylist]
  best_list = []
  best_dist = []
  for i in range(0,len(mylist)):
    counts[i] == maxcount
    best_list.append(mylist[i])
    best_dist.append(dist[i])
  abs_best_dist = min(best_dist)
  for i in range(0,len(best_dist)):
    if best_dist[i] == abs_best_dist: 
      return best_list[i]
  sys.stderr.write("Warning: trouble finding best\n")
  return best_list[0]

class Transcriptome:
  """a class to store a transcriptome

  :param gpd_file: filename
  :param ref_fasta:
  :type gpd_file: string
  :type ref_fasta: dict()
  """
  def __init__(self,gpd_file=None,ref_fasta=None):
    self.transcripts = []
    if gpd_file:
      from seqtools.format.GPD import GPD
      with open(gpd_file) as inf:
        for line in inf:
          self.transcripts.append(GPD(line))
    if ref_fasta:
      for i in range(0,len(self.transcripts)):
        self.transcripts[i].get_sequence(ref_fasta)
  def dump_serialized(self):
    sx = base64.b64encode(zlib.compress(pickle.dumps([x.dump_serialized() for x in self.transcripts])))
    return sx
  def load_serialized(self,instr):
    txs = []
    for v in pickle.loads(zlib.decompress(base64.b64decode(instr))):
      tx = Transcript()
      tx.load_serialized(v)
      txs.append(tx)
    self.transcripts = txs

  def get_transcripts(self):
    return self.transcripts
      
  def add_transcript(self,transcript):
    self.transcripts.append(transcript)

  def __str__(self):
    ostr = ''
    ostr += "Transcriptome containing "+str(len(self.transcripts))+" transcripts "
    ostr += "covering "+str(sum([x.get_length() for x in self.transcripts]))+" bases"
    return ostr
