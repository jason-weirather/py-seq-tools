import sys
import seqtools.structure.transcript
from seqtools.graph import Graph, Node, Edge
from seqtools.range import GenomicRange
from seqtools.statistics import mode
from random import shuffle

class ExonGroup:
   """ Class for accessing a group of exons that likely represent the same exon"""
   def __init__(self):
      self._exons = []
   def add_exon(self,ex):
      #for preex in self._exons:
      #   if not preex.overlaps(ex): raise ValueError('non overlapping exon')
      self._exons.append(ex)

   @property
   def chr(self): 
      if len(self._exon)==0: return None
      return self._exons[0].chr

   @property
   def start_range(self):
      """Similar to the junction range but don't need to check for leftmost or rightmost"""
      if len(self._exons) == 0: return None
      return GenomicRange(self._exons[0].chr,
             min([x.start for x in self._exons]),# must be part of junction
             max([x.start for x in self._exons]))
   @property
   def end_range(self):
      """Similar to the junction range but don't need to check for leftmost or rightmost"""
      if len(self._exons) == 0: return None
      return GenomicRange(self._exons[0].chr,
             min([x.end for x in self._exons]),
             max([x.end for x in self._exons]))

   @property
   def left_junction_range(self):
      if len(self._exons) == 0: return None
      return GenomicRange(self._exons[0].chr,
             min([x.end for x in self._exons if not x.rightmost]),# must be part of junction
             max([x.end for x in self._exons if not x.rightmost]))
   @property
   def right_junction_range(self):
      if len(self._exons) == 0: return None
      return GenomicRange(self._exons[0].chr,
             min([x.start for x in self._exons if not x.leftmost]),
             max([x.start for x in self._exons if not x.leftmost]))

   def consensus(self,position=None):
      if position == 'leftmost':
         return GenomicRange(self._exons[0].chr,
            min([x.start for x in self._exons]),
            mode([x.end for x in self._exons if not x.rightmost]))
      if position == 'rightmost':
         return GenomicRange(self._exons[0].chr,
            mode([x.start for x in self._exons if not x.leftmost]),
            max([x.end for x in self._exons]))
      if position == 'internal':
         return GenomicRange(self._exons[0].chr,
            mode([x.start for x in self._exons if not x.leftmost]),
            mode([x.end for x in self._exons if not x.rightmost]))
      if position == 'single':
         return GenomicRange(self._exons[0].chr,
            min([x.start for x in self._exons]),
            max([x.end for x in self._exons]))
      return GenomicRange(self._exons[0].chr,
            mode([x.start for x in self._exons]),
            mode([x.end for x in self._exons]))


class FuzzyTranscript(seqtools.structure.transcript.Transcript):
   """A transcript that has one single set of exons, but the bounds may differ"""
   def __init__(self,initial_transcripts,tolerance=0,evidence=1):
      self._evidence = evidence
      super(FuzzyTranscript,self).__init__([])
      if not self.direction: self.set_strand(initial_transcripts[0].strand)
      """Initialize with one or more transcripts of the same layout"""
      self._initial = initial_transcripts
      self._transcripts = self._initial[:]
      self._tolerance = tolerance
      """Sanity check the initial for compatibility"""
      self._num = initial_transcripts[0].get_exon_count()
      self._exon_groups = [ExonGroup() for x in range(0,self._num)]
      for tx in self._initial:
         for i in range(0,self._num):
            self._exon_groups[i].add_exon(tx.exons[i])
      """Now exon groups is set up"""

   #@property
   #def direction(self):
   #   return self._transcripts[0].direction
   #@property
   #def strand(self):  return self.direction

   @property
   def _rngs(self):
      """This is where we should also enforce evidence requirements"""
      outputs = []
      if len(self._exon_groups)==1:
         return [self._exon_groups.consensus('single')]
      z = 0 #output count
      begin = 0
      meeting_criteria = [i for i in range(0,len(self._exon_groups)) if  len(self._exon_groups) >= self._evidence]
      if len(meeting_criteria) == 0: return []
      finish = len(meeting_criteria)
      if len(meeting_criteria) > 0:
         begin = meeting_criteria[0]
         finish = meeting_criteria[-1]
      for i in range(0,len(self._exon_groups)):
         if z == begin:
            outputs.append(self._exon_groups[i].consensus('leftmost'))
         elif z == finish:
            #len(self._exon_groups)-1:
            outputs.append(self._exon_groups[i].consensus('rightmost'))
         else:
            outputs.append(self._exon_groups[i].consensus('internal'))
         z += 1
      v = [seqtools.structure.transcript.Exon(x) for x in outputs]
      v[0].set_leftmost()
      v[-1].set_rightmost()
      return v


   @property
   def exons(self):
      return self._rngs

   def add_transcripts(self,txs):
      """We traverse through the other transcripts and try to add to these groups
      """
      passed = []
      for tx2 in txs:
         for tx1 in self._initial:
            jov = tx1.junction_overlap(tx2,self._tolerance)
            sub = jov.is_subset()
            if sub == 1 or sub == 2:
               passed.append(tx2)
               break
      if len(passed) == 0:
         sys.stderr.write("Warning unable to add\n")
         return
      for tx in txs:
         self.add_transcript(tx)
      return

   def add_transcript(self,tx):
      """add a single transcript"""
      candidates = tx.junctions
      targets = self.junctions
      matches = []
      for i in range(0,len(targets)):
         for j in range(0,len(candidates)):
            if targets[i].overlaps(candidates[j],self._tolerance):
               matches.append([i,j])
      if len(matches) != len(candidates): return
      if len(matches) > 1:
         if False in  [(matches[i+1][0]-matches[i][0])==1 and
                       (matches[i+1][1]-matches[i][1])==1 for i in range(0,len(matches)-1)]:
            return
      # nowe we can add them
      for m in matches:
         self._exon_groups[m[0]].add_exon(tx.exons[m[1]])
         self._exon_groups[m[0]+1].add_exon(tx.exons[m[1]+1])
      self._transcripts.append(tx)

   @property
   def transcripts(self):
      return self._transcripts

class CompatibleGraph(FuzzyTranscript):
   """A group of transcripts representing a single transcript
      
      Its like a fuzzy transcript, except its created from a graph
   """
   def __init__(self,graph,tolerance=0,downsample=None,evidence=1):
      self._graph = graph
      self._tolerance = tolerance
      """One sanity check is to make sure each code contains the same junction count"""
      for x in [n.payload_list for n in self._graph.nodes]:
         ecnts = [y.get_exon_count() for y in x]
         if min(ecnts) != max(ecnts):
            raise ValueError('compatible nodes should not be having different junction counts')
      """traverse the graph and build a transcript"""
      root = self._graph.roots[0] #there should only be one root of a directed graph
      children = self._graph.get_children(root)
      """Start at the root transcripts"""
      initial_transcripts = root.payload_list
      """Now we can initialize as a fuzzy transcript"""
      super(CompatibleGraph,self).__init__(initial_transcripts,tolerance,evidence)
      #ft = FuzzyTranscript(initial_transcripts,self._tolerance)
      all_transcripts = []
      all_transcripts += [x.payload_list for x in children]
      if downsample:
         all_transcripts = do_downsample(all_transcripts,downsample)
      for transcripts in all_transcripts:
         self.add_transcripts(transcripts)

def do_downsample(txs,cnt):
   v = txs[:]
   shuffle(v)
   return v[0:cnt]


class Deconvolution:
   """Take a group of transcripts.  
      Preferably just a locus worth and separate them 
      into transcripts representative of each group.
   """
   def __init__(self,txs=[]):
      self._transcripts = txs
   def add_transcript(self,tx):
      self._transcripts.append(tx)
   def parse(self,tolerance=0,downsample=None,evidence=2):
      """Divide out the transcripts.  allow junction tolerance if wanted"""
      g = Graph()
      nodes = [Node(x) for x in self._transcripts]
      for n in nodes: g.add_node(n)
      for i in range(0,len(nodes)):
         for j in range(0,len(nodes)):
            if i == j: continue
            jov = nodes[i].payload.junction_overlap(nodes[j].payload,tolerance)
            sub = jov.is_subset()
            if not sub: continue
            if sub == 1:
               g.add_edge(Edge(nodes[i],nodes[j]))
               g.add_edge(Edge(nodes[j],nodes[i]))
            if sub == 2:
               g.add_edge(Edge(nodes[i],nodes[j]))
      g.merge_cycles()
      roots = g.roots
      groups = []
      for r in roots:
         g2 = g.get_root_graph(r)
         c = CompatibleGraph(g2,tolerance,downsample,evidence)
         groups.append(c)
      return groups
