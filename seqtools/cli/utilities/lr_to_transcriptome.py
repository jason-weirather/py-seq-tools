"""CLI script to annotate gpd entries
"""

import sys, argparse, gzip, re, itertools, os, uuid
from multiprocessing import cpu_count
from tempfile import mkdtemp, gettempdir
from random import shuffle
from seqtools.format.gpd import GPD, GPDStream
from seqtools.structure.transcriptome import Transcriptome
from seqtools.structure.transcript import Transcript
from seqtools.stream import OrderedStream
from seqtools.stream import MultiLocusStream
from seqtools.structure.transcript.group import Deconvolution
from seqtools.graph import Graph, Node, Edge
from seqtools.range import GenomicRange
from subprocess import Popen, PIPE


def main(args):
   """Read any reference transcriptome we can"""
   txome = Transcriptome()
   #read the reference gpd if one is gven
   if args.reference: 
      rinf = None
      if re.search('\.gz$',args.reference):
         rinf = gzip.open(args.reference)
      else:
         rinf = open(args.reference)
      sys.stderr.write("Reading in reference\n")
      z = 0
      # populate txome with reference transcripts for each chromosome
      for line in rinf:
         z += 1
         gpd = GPD(line)
         gpd.set_payload(z)
         if z%100 == 0:  sys.stderr.write(str(z)+"          \r")
         txome.add_transcript(gpd)
      rinf.close()
      sys.stderr.write(str(z)+"          \r")
      sys.stderr.write("\n")
   txome.sort_transcripts()
   sys.stderr.write("Buffering mappings\n")
   inf = sys.stdin
   if args.input != '-':
      if re.search('\.gz$',args.input):
         inf = gzip.open(args.input)
      else:
         inf = open(args.input)
   tof = gzip.open(args.tempdir+'/reads.gpd.gz','w')
   for line in inf: tof.write(line.rstrip()+"\n")
   tof.close()
  
   sys.stderr.write("1. Process by overlapping locus.\n")
   annotated_reads(txome,gzip.open(args.tempdir+'/reads.gpd.gz'),'initial',args)

   #of.close()

   """Now sequences we should be able to annotate as partial have
      consensus transcripts.  Lets put those to annotation scrutiny again.
   """
   sys.stderr.write("2. Reannotated partially annotated by overlapping locus.\n")
   annotated_reads(txome,gzip.open(args.tempdir+'/initial/partial_annotated_multiexon.sorted.gpd.gz'),'repartial',args)


   txs = {}
   for transcript in txome.transcripts:
      txs[transcript.name] = transcript
   detected = {}
   tinf = gzip.open(args.tempdir+'/initial/annotated.txt.gz')
   for line in tinf:
      f = line.rstrip().split("\t")
      detected[f[1]] = True
   tinf.close()
   tinf = gzip.open(args.tempdir+'/repartial/annotated.txt.gz')
   remove_partial = {}
   for line in tinf:
      f = line.rstrip().split("\t")
      detected[f[1]] = True
      remove_partial[f[1]] = True
   tinf.close()
   tof = gzip.open(args.tempdir+'/candidate.gpd.gz','w')
   for name in detected:
      tof.write(txs[name].get_gpd_line()+"\n")
   tinf = gzip.open(args.tempdir+'/initial/partial_annotated_multiexon.sorted.gpd.gz')
   for line in tinf:
      f = line.rstrip().split("\t")
      if f[1] not in remove_partial:
         tof.write(line)
   tinf = gzip.open(args.tempdir+'/initial/unannotated_singleexon.sorted.gpd.gz')
   for line in tinf:
      f = line.rstrip().split("\t")
      f[0] = str(uuid.uuid4())
      tof.write("\t".join(f)+"\n")
   tinf.close()
   tof.close()
   sort_gpd(args.tempdir+'/candidate.gpd.gz',args.tempdir+'/candidate.sorted.gpd.gz',args)
   """We can ignore the partial annotations that have been detected"""
   ntxome = Transcriptome()
   tinf = gzip.open(args.tempdir+'/candidate.sorted.gpd.gz')
   for line in tinf: ntxome.add_transcript(GPD(line))     
   annotated_reads(ntxome,gzip.open(args.tempdir+'/initial/unannotated_multiexon.sorted.gpd.gz'),'unannotated',args)
   """now we know which unannotated reads actually have annotations"""
   tinf.close()

   tof = gzip.open(args.tempdir+'/final.gpd.gz','w')
   tinf = gzip.open(args.tempdir+'/candidate.sorted.gpd.gz')   
   for line in tinf: tof.write(line)
   tinf.close()
   tinf = gzip.open(args.tempdir+'/unannotated/unannotated_multiexon.sorted.gpd.gz')   
   for line in tinf: tof.write(line)
   tinf.close()
   tinf = gzip.open(args.tempdir+'/unannotated/partial_annotated_multiexon.sorted.gpd.gz')   
   for line in tinf: tof.write(line)
   tinf.close()
   tof.close()
   sort_gpd(args.tempdir+'/final.gpd.gz',args.tempdir+'/final.sorted.gpd.gz',args)
   """Prepare outputs"""
   of = sys.stdout
   if args.output:
      if re.search('\.gz$',args.output):
         of = gzip.open(args.output,'w')
      else:
         of = open(args.output,'w')
   tinf = gzip.open(args.tempdir+'/final.sorted.gpd.gz')
   for line in tinf: of.write(line)
   tinf.close()
   of.close()   
      

def annotated_reads(txome,sorted_read_stream,output,args):
   if not os.path.exists(args.tempdir+'/'+output):
      os.makedirs(args.tempdir+'/'+output)
   ts = OrderedStream(iter(txome.transcript_stream()))
   rs = OrderedStream(GPDStream(sorted_read_stream))
   mls = MultiLocusStream([ts,rs])
   
   aof = gzip.open(args.tempdir+'/'+output+'/annotated.txt.gz','w')
   of1 = gzip.open(args.tempdir+'/'+output+'/unannotated_multiexon.gpd.gz','w')
   of2 = gzip.open(args.tempdir+'/'+output+'/partial_annotated_multiexon.gpd.gz','w')
   of3 = gzip.open(args.tempdir+'/'+output+'/unannotated_singleexon.gpd.gz','w')
   z = 0
   for ml in mls:
      refs, reads = ml.payload
      if len(refs) == 0: continue
      """Check and see if we have single exons annotated"""
      single_exon_refs = [x for x in refs if x.get_exon_count()==1]
      single_exon_reads = [x for x in reads if x.get_exon_count()==1]
      #print str(len(single_exon_refs))+" "+str(len(single_exon_reads))
      unannotated_single_exon_reads = []
      for seread in single_exon_reads:
         ovs  = [(x.exons[0].length,seread.exons[0].length,x.exons[0].overlap_size(seread.exons[0]),x) for x in single_exon_refs if x.exons[0].overlaps(seread.exons[0])]
         #check for minimum overlap
         ovs = [x for x in ovs if x[2] >= args.single_exon_minimum_overlap]
         ovs = [x for x in ovs if float(x[2])/float(max(x[0],x[1])) >= args.single_exon_mutual_overlap]
         ovs = sorted(ovs,key=lambda x: float(x[2])/float(max(x[0],x[1])),reverse=True)
         if len(ovs) == 0: 
            unannotated_single_exon_reads.append(seread)
            continue
         best_ref = ovs[0][3]
         aof.write(seread.name+"\t"+best_ref.name+"\tSE"+"\n")

      """May want an optional check for any better matches among exons"""
      #now we can look for best matches among multi-exon transcripts
      reads = [x for x in reads if x.get_exon_count() > 1]
      multiexon_refs = [x for x in refs if x.get_exon_count() > 1]
      unannotated_multi_exon_reads = []
      partial_annotated_multi_exon_reads = []
      for read in reads:
         # we dont' need to have multiple exons matched. one is enough to call
         ovs = [y for y in [(x,x.exon_overlap(read,
                                              multi_minover=args.multi_exon_minimum_overlap,
                                              multi_endfrac=args.multi_exon_end_frac,
                                              multi_midfrac=args.multi_exon_mid_frac,
                                              multi_consec=False)) for x in multiexon_refs] if y[1]]
         for o in ovs: o[1].analyze_overs()
         full =  sorted([x for x in ovs if x[1].is_subset()==1],
                        key = lambda y: float(y[1].overlap_size())/float(max(y[1].tx_obj1.length,y[1].tx_obj2.length)),
                        reverse=True
                       )
         if len(full) > 0:
            aof.write(read.name+"\t"+full[0][0].name+"\tfull"+"\n")
            continue
         subset =  sorted([x for x in ovs if x[1].is_subset()==2],
                          key = lambda y: (y[1].match_exon_count(),
                                           y[1].min_overlap_fraction()),
                          reverse = True
                         )
         if len(subset) > 0:
            aof.write(read.name+"\t"+subset[0][0].name+"\tpartial"+"\n")
            continue
         #check for supersets
         superset = sorted([x for x in ovs if x[1].is_subset()==3],
                          key = lambda y: (y[1].match_exon_count(),
                                           y[1].min_overlap_fraction()),
                          reverse = True
                         )
         if len(superset) > 0:
            partial_annotated_multi_exon_reads.append((read,superset[0][0]))
            #print read.name+"\t"+superset[0][0].name+"\tsuper"
            continue
         #check for noncompatible overlaps
         overset = sorted([x for x in ovs if x[1].match_exon_count > 0],
                          key = lambda y: (y[1].consecutive_exon_count(),
                                           y[1].min_overlap_fraction()),
                          reverse = True
                         )
         #print [(x[1].consecutive_exon_count(), x[1].min_overlap_fraction()) for x in overset]
         if len(overset) > 0:
            partial_annotated_multi_exon_reads.append((read,overset[0][0]))
            #print read.name+"\t"+overset[0][0].name+"\tover"
            continue
         unannotated_multi_exon_reads.append(read)
      """Now we have partially annotated multi and unannotated multi and unannotated single"""
      if len(unannotated_multi_exon_reads) > 0:
         sys.stderr.write(str(z)+" "+str(len(unannotated_multi_exon_reads))+"   \r")
         d = Deconvolution(downsample(unannotated_multi_exon_reads,args.downsample_locus))
         groups = d.parse(tolerance=20,downsample=args.downsample)
         for tx in groups:
            z+=1
            of1.write(tx.get_gpd_line()+"\n")
      if len(partial_annotated_multi_exon_reads) > 0:
         sys.stderr.write(str(z)+" "+str(len(partial_annotated_multi_exon_reads))+"   \r")
         ### set the direction of the transcript
         for v in partial_annotated_multi_exon_reads:
            v[0].set_strand(v[1].direction)
         ### set the gene name of the transcript
         for v in partial_annotated_multi_exon_reads:
            v[0].set_gene_name(v[1].gene_name)         
         d = Deconvolution(downsample([x[0] for x in partial_annotated_multi_exon_reads],args.downsample_locus))
         groups = d.parse(tolerance=20,downsample=args.downsample,use_gene_names=True)
         for tx in groups:
            z += 1
            of2.write(tx.get_gpd_line()+"\n")
      """Do the unannotated single exon reads"""
      g = Graph()
      for r in unannotated_single_exon_reads:
         if len([x for x in partial_annotated_multi_exon_reads if x[0].overlaps(r)]) > 0: continue
         if len([x for x in unannotated_multi_exon_reads if x.overlaps(r)]) > 0: continue
         if len([x for x in txome.transcripts if x.overlaps(r)]) > 0: continue
         g.add_node(Node(r))
      for i in range(0,len(g.nodes)):
         for j in range(0,len(g.nodes)):
            if i == j: continue
            if g.nodes[i].payload.overlaps(g.nodes[j].payload):
               g.add_edge(Edge(g.nodes[i],g.nodes[j]))
      g.merge_cycles()
      for r in g.roots:
         se = []
         se += r.payload_list
         children = g.get_children(r)
         for child in children:  se += child.payload_list
         rng = GenomicRange(se[0].exons[0].chr,
                            min([x.exons[0].start for x in se]),
                            max([x.exons[0].end for x in se]))
         tx = Transcript([rng],Transcript.Options(direction='+'))
         of3.write(tx.get_gpd_line()+"\n")

   of1.close()
   of2.close()
   of3.close()
   sys.stderr.write("\n")      
   """Sort our progress"""
   sys.stderr.write("sort transcriptome\n")
   sort_gpd(args.tempdir+'/'+output+'/partial_annotated_multiexon.gpd.gz',args.tempdir+'/'+output+'/partial_annotated_multiexon.sorted.gpd.gz',args)
   sort_gpd(args.tempdir+'/'+output+'/unannotated_multiexon.gpd.gz',args.tempdir+'/'+output+'/unannotated_multiexon.sorted.gpd.gz',args)
   sort_gpd(args.tempdir+'/'+output+'/unannotated_singleexon.gpd.gz',args.tempdir+'/'+output+'/unannotated_singleexon.sorted.gpd.gz',args)
   """We still have the unannotated single exon reads to deal with"""


def sort_gpd(infile,outfile,args):
   tinf = gzip.open(infile)
   cmd1 = 'seq-tools sort --gpd - '
   cmd2 = 'gzip'
   if args.threads > 1: cmd1 += ' --threads '+str(args.threads)
   tof = open(outfile,'w')
   p2 = Popen(cmd2.split(),stdin=PIPE,stdout=tof)
   p1 = Popen(cmd1.split(),stdin=PIPE,stdout=p2.stdin)
   for line in tinf: p1.stdin.write(line)
   p1.communicate()
   p2.communicate()
   tinf.close()
   tof.close()   


def downsample(txs,cnt):
   if cnt >= len(txs): return txs
   v = txs[:]
   shuffle(v)
   return v[0:cnt]

def setup_tempdir(args):
  if args.specific_tempdir:
    if not os.path.exists(args.specific_tempdir):
      os.makedirs(args.specific_tempdir.rstrip('/'))
    args.tempdir = args.specific_tempdir.rstrip('/')
    if not os.path.exists(args.specific_tempdir.rstrip('/')):
      sys.stderr.write("ERROR: Problem creating temporary directory\n")
      sys.exit()
  else:
    args.tempdir = mkdtemp(prefix="weirathe.",dir=args.tempdir.rstrip('/'))
    if not os.path.exists(args.tempdir.rstrip('/')):
      sys.stderr.write("ERROR: Problem creating temporary directory\n")
      sys.exit()
  if not os.path.exists(args.tempdir):
    sys.stderr.write("ERROR: Problem creating temporary directory\n")
    sys.exit()
  return



def do_inputs():
  parser = argparse.ArgumentParser(description="build a transcriptome from long reads",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="Use - for STDIN ordered GPD")
  parser.add_argument('-o','--output',help="output file otherwise STDOUT")
  parser.add_argument('--threads',type=int,default=cpu_count(),help="Number of threads to convert names")
  parser.add_argument('-r','--reference',help="reference gpd")
  parser.add_argument('--downsample',type=int,default=50,help='at most this many reads')
  parser.add_argument('--downsample_locus',type=int,default=50,help="limit the number of locus transcripts to this amount")
  group1 = parser.add_argument_group('Single Exon Parameters')
  group1.add_argument('--single_exon_minimum_overlap',type=int,default=20,help="minimum bases a read must overlap with a reference single exon transcript")
  group1.add_argument('--single_exon_mutual_overlap',type=int,default=0.2,help="minimum bases a read must overlap with a reference single exon transcript")
  group2 = parser.add_argument_group('Multi Exon Parameters')
  group2.add_argument('--multi_exon_minimum_overlap',type=int,default=10,help="minimum amount exons must overlap")
  group2.add_argument('--multi_exon_end_frac',type=float,default=0,help="minimum mutual overlap of end exons")
  group2.add_argument('--multi_exon_mid_frac',type=float,default=0.8,help="minimum mutual overlap of middle exons")
  group = parser.add_mutually_exclusive_group()
  group.add_argument('--tempdir',default=gettempdir(),help="The temporary directory is made and destroyed here.")
  group.add_argument('--specific_tempdir',help="This temporary directory will be used, but will remain after executing.")
  args = parser.parse_args()
  setup_tempdir(args)
  return args
  
def external_cmd(cmd):
  cache_argv = sys.argv
  sys.argv = cmd
  args = do_inputs()
  main(args)
  sys.argv = cache_argv

if __name__=="__main__":
  args = do_inputs()
  main(args)
