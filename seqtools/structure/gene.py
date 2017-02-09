"""a module for working with families of transcripts

"""
import sys, random, string, uuid
from collections import namedtuple
from seqtools.range import GenomicRange
from seqtools.range.multi import ranges_to_coverage, merge_ranges
from seqtools.sequence import rc
import seqtools.graph
from math import sqrt


class TranscriptLoci:
  """combine together compatible multiple transcript groups to form
   a simpler set of transcripts """
  def __init__(self):
    #self.transcripts = []
    self.merge_rules = TranscriptLociMergeRules('is_any_overlap')
    self.merge_rules.set_juntol(10)
    self.g = seqtools.graph.Graph()   

  def __str__(self):
    return str(len(self.g.get_nodes()))+ " nodes"  

  def remove_transcript(self,tx_id):
    """Remove a transcript from the locus by its id

    :param tx_id:
    :type tx_id: string
    """
    txs = self.get_transcripts()
    if tx_id not in [x.id for x in txs]:
      return
    tx = [x for x in txs if x.id==tx_id][0]
    for n in [x for x in self.g.get_nodes()]:
      if tx_id not in [y.id for y in n.payload]:
        continue
      n.payload.remove(tx)
      if len(n.payload)==0:
        self.g.remove_node(n)      
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
      for ex in [x.range for x in tx.exons]: bedarray.append(ex)
    cov = ranges_to_coverage(bedarray)
    results = {}
    for tx in self.get_transcripts():
      tlen = tx.range.length
      bcov = []
      for ex in [x.range for x in tx.exons]:     
        excov = [[x.overlap_size(ex),x.payload] for x in cov]
        for coved in [x for x in excov if x[0] > 0]:
          bcov.append(coved)
      total_base_coverage = sum([x[0]*x[1] for x in bcov])
      average_coverage = float(total_base_coverage)/float(tlen)
      minimum_bases_covered = sum([x[0] for x in bcov if x[1] >= mindepth])
      fraction_covered_at_minimum = float(minimum_bases_covered)/float(tlen)
      res = {'tx':tx,'average_coverage':average_coverage,'fraction_covered':fraction_covered_at_minimum,'mindepth':mindepth,'length_covered':minimum_bases_covered}
      results[tx.id] = res
      #print average_coverage
      #print fraction_covered_at_minimum
      #print tlen
      #tcov = float(bcov)/float(tlen)
      #print tcov
    #for c in cov:
    #  print c
    return results

  @property
  def range(self):
    """Return the range the transcript loci covers

    :return: range
    :rtype: GenomicRange
    """
    chrs = set([x.range.chr for x in self.get_transcripts()])
    if len(chrs) != 1: return None
    start = min([x.range.start for x in self.get_transcripts()])
    end = max([x.range.end for x in self.get_transcripts()])
    return GenomicRange(list(chrs)[0],start,end)

  def get_transcripts(self):
    """ a list of the transcripts in the locus"""
    txs = []
    for pays in [x.payload for x in self.g.get_nodes()]:
      for pay in pays:
        txs.append(pay)
    return txs

  def partition_loci(self,verbose=False):
    """ break the locus up into unconnected loci

    :return: list of loci
    :rtype: TranscriptLoci[]
    """
    self.g.merge_cycles()
    #sys.stderr.write(self.g.get_report()+"\n")
    gs = self.g.partition_graph(verbose=verbose)
    tls = [] # makea list of transcript loci
    for g in gs:
      tl = TranscriptLoci()
      tl.merge_rules = self.merge_rules
      ns = g.get_nodes()
      for n in [x.payload for x in ns]:
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
    for y in [x.payload for x in self.g.get_nodes()]:
      if tx.id in [z.id for z in y]:
        sys.stderr.write("WARNING tx is already in graph\n")
        return True
    # transcript isn't part of graph yet
    n = seqtools.graph.Node([tx])

    other_nodes = self.g.get_nodes()
    self.g.add_node(n)
    # now we need to see if its connected anywhere
    for n2 in other_nodes:
     tx2s = n2.payload
     for tx2 in tx2s:
      # do exon overlap
      er = self.merge_rules.get_exon_rules()
      # if we are doing things by exon
      if (self.merge_rules.get_use_single_exons() and (tx.get_exon_count() == 1 or tx2.get_exon_count() == 1)) or \
         (self.merge_rules.get_use_multi_exons() and (tx.get_exon_count() > 1 and tx2.get_exon_count() > 1)):
        eo = tx.exon_overlap(tx2,multi_minover=er['multi_minover'],multi_endfrac=er['multi_endfrac'],multi_midfrac=er['multi_midfrac'],single_minover=er['single_minover'],single_frac=er['single_frac'])
        if self.merge_rules.get_merge_type() == 'is_compatible':
          if eo.is_compatible():
            self.g.add_edge(seqtools.graph.Edge(n,n2),verbose=False)
            self.g.add_edge(seqtools.graph.Edge(n2,n),verbose=False)
        elif self.merge_rules.get_merge_type() == 'is_subset':
          r = eo.is_subset()
          if r == 2 or r == 1:
            self.g.add_edge(seqtools.graph.Edge(n,n2),verbose=False)
          if r == 3 or r == 1:
            self.g.add_edge(seqtools.graph.Edge(n2,n),verbose=False)
        elif self.merge_rules.get_merge_type() == 'is_full_overlap':
          if eo.is_full_overlap():
            self.g.add_edge(seqtools.graph.Edge(n,n2),verbose=False)
            self.g.add_edge(seqtools.graph.Edge(n2,n),verbose=False)
        elif self.merge_rules.get_merge_type() == 'is_any_overlap':
          if eo.match_exon_count() > 0:
            self.g.add_edge(seqtools.graph.Edge(n,n2),verbose=False)
            self.g.add_edge(seqtools.graph.Edge(n2,n),verbose=False)        
            
      if self.merge_rules.get_use_junctions():
        # do junction overlap
        jo = tx.junction_overlap(tx2,self.merge_rules.get_juntol())
        #print jo.match_junction_count()
        if self.merge_rules.get_merge_type() == 'is_compatible':
          if jo.is_compatible():
            self.g.add_edge(seqtools.graph.Edge(n,n2),verbose=False)
            self.g.add_edge(seqtools.graph.Edge(n2,n),verbose=False)
        elif self.merge_rules.get_merge_type() == 'is_subset':
          r = jo.is_subset()
          if r == 2 or r == 1:
            self.g.add_edge(seqtools.graph.Edge(n,n2),verbose=False)
          if r == 3 or r == 1:
            self.g.add_edge(Seqtools.graph.Edge(n2,n),verbose=False)
        elif self.merge_rules.get_merge_type() == 'is_full_overlap':
          if jo.is_full_overlap():
            self.g.add_edge(seqtools.graph.Edge(n,n2),verbose=False)
            self.g.add_edge(seqtools.graph.Edge(n2,n),verbose=False)
        elif self.merge_rules.get_merge_type() == 'is_any_overlap':
          if jo.match_junction_count() > 0:
            self.g.add_edge(seqtools.graph.Edge(n,n2),verbose=False)
            self.g.add_edge(seqtools.graph.Edge(n2,n),verbose=False)        
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
      leftcoord = min([x.exons[0].range.start for x in self.transcripts])
      rightcoord = max([x.exons[-1].range.end for x in self.transcripts])
      e = Exon(GenomicRange(x.exons[0].range.chr,leftcoord,rightcoord))
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
                               min([x.range.start for x in left_exons]),
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
                               max([x.range.end for x in right_exons])))
    e_right.set_left_junc(out.junctions[-1])
    out.exons.append(e_right)
    return out

  def add_transcript(self,tx,juntol=0,verbose=True):
    if tx.id in self._transcript_ids: return True
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
    self._transcript_ids.add(tx.id)
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
    ostr += "covering "+str(sum([x.range.length for x in self.transcripts]))+" bases"
    return ostr

def trim_ordered_range_list(ranges,start,finish):
  """A function to help with slicing a mapping
     Start with a list of ranges and get another list of ranges constrained by start (0-indexed) and finish (1-indexed)

     :param ranges: ordered non-overlapping ranges on the same chromosome
     :param start: start 0-indexed
     :param finish: ending 1-indexed
     :type ranges: GenomicRange []
     :type start: Int
     :type finish: Int
     :return: non-overlapping ranges on same chromosome constrained by start and finish
     :rtype: GenomicRange []
  """
  z = 0
  keep_ranges = []
  for inrng in self.ranges:
    z+=1
    original_rng = inrng
    rng = inrng.copy() # we will be passing it along and possibly be cutting it
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
  return keep_ranges
