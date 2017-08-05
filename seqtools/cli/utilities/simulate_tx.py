"""CLI script to simulate transcriptome reads
"""

import sys, argparse, gzip, os, itertools
from seqtools.format.gpd import GPDStream
from seqtools.format.fasta import FASTAData
from seqtools.structure.transcriptome import Transcriptome
from seqtools.simulation.emitter import TranscriptomeEmitter
from seqtools.simulation.emitter import ReadEmitter
from seqtools.simulation.randomsource import RandomSource
from multiprocessing import cpu_count, Pool
from seqtools.simulation.permute import ErrorMakerFlatRate
from seqtools.quantification import TPMCalculator

def main(args):
   rand = None
   if args.seed:
      rand = RandomSource(args.seed)
   else:
      rand = RandomSource()
   if args.transcriptome[-3:] == '.gz':
      gpd = GPDStream(gzip.open(args.transcriptome))
   else:
      gpd = GPDStream(open(args.transcriptome))
   if args.genome[-3:] == '.gz':
      fasta = FASTAData(gzip.open(args.genome).read())
   else:
      fasta = FASTAData(open(args.genome).read())
   txome = Transcriptome(gpd,fasta)
   if not os.path.exists(args.output):
      os.makedirs(args.output)

   if args.short_read_count:
      gen = prepare_gen(args,txome,rand,20000,args.short_read_count)
      calc = TPMCalculator(txome)
      if args.threads > 1: p = Pool(processes=args.threads)
      of1 = gzip.open(args.output+'/SR_output_1.fq.gz','w')
      of2 = gzip.open(args.output+'/SR_output_2.fq.gz','w')
      of3 = gzip.open(args.output+'/SR_original_1.fa.gz','w')
      of4 = gzip.open(args.output+'/SR_original_2.fa.gz','w')
      of5 = gzip.open(args.output+'/SR_original_1.gpd.gz','w')
      of6 = gzip.open(args.output+'/SR_original_2.gpd.gz','w')
      z = 0
      #gen = prepare_gen(args,gpd,fasta,rand)
      if args.threads > 1:
         func = p.imap(get_short_reads,gen)
      else:
         func = itertools.imap(get_short_reads,gen)
      for reads in func: 
         """Emit the short reads first"""
         for e in reads:
            z += 1
            if z%100==0: sys.stderr.write(str(z)+"     \r")
            if z > args.short_read_count: break
            of1.write(str(e.reads.left))
            of2.write(str(e.reads.right))
            of3.write(">"+e.reads.left.name+"\tleft\t"+e.transcript.name+"\t"+str(e.source.fragment.length)+"\n"+str(e.source.left.sequence)+"\n")
            of4.write(">"+e.reads.right.name+"\tright\t"+e.transcript.name+"\t"+str(e.source.fragment.length)+"\n"+str(e.source.right.sequence)+"\n")
            of5.write(e.source.left.get_gpd_line(transcript_name=e.reads.left.name,gene_name=e.reads.left.name)+"\n")
            of6.write(e.source.right.get_gpd_line(transcript_name=e.reads.right.name,gene_name=e.reads.right.name)+"\n")
            # save the tpm value
            calc.add_count(e.transcript.name)
      sys.stderr.write("\n")
      of1.close()
      of2.close()
      of3.close()
      of4.close()
      of5.close()
      of6.close()
      sys.stderr.write("\n")
      of7 = gzip.open(args.output+'/SR_output.quant.gz','w')
      of7.write("name\tTPM\tFPKM\tcount\n")
      for name in sorted([y.name for y in txome.transcripts],key=lambda x: calc.TPM(x), reverse=True):
         of7.write(name+"\t"+str(calc.TPM(name))+"\t"+str(calc.FPKM(name))+"\t"+str(calc.count(name))+"\n")
      of7.close()
      #if args.threads > 1: 
      #   p.close()
      #   p.join()

   if args.long_read_count:
      gen = prepare_gen(args,txome,rand,1000,args.long_read_count)
      calc = TPMCalculator(txome)
      if args.threads > 1: p = Pool(processes=args.threads)
      of1 = gzip.open(args.output+'/LR_output.fq.gz','w')
      of2 = gzip.open(args.output+'/LR_original.fa.gz','w')
      of3 = gzip.open(args.output+'/LR_original.gpd.gz','w')
      z = 0
      if args.threads > 1:
         func = p.imap(get_long_reads,gen)
      else:
         func = itertools.imap(get_long_reads,gen)
      for reads in func: 
         for e in reads:
            z += 1
            if z%100==0: sys.stderr.write(str(z)+"     \r")
            if z > args.long_read_count: break
            of1.write(str(e.reads.fragment))
            of2.write(">"+e.reads.fragment.name+"\tlong\t"+e.transcript.name+"\t"+str(e.source.fragment.length)+"\n"+str(e.source.fragment.sequence)+"\n")
            of3.write(e.source.fragment.get_gpd_line(transcript_name=e.reads.fragment.name,gene_name=e.reads.fragment.name)+"\n")
            calc.add_count(e.transcript.name)
      sys.stderr.write("\n")
      of1.close()
      of2.close()
      of3.close()
      of4 = gzip.open(args.output+'/LR_output.quant.gz','w')
      of4.write("name\tTPM\tFPKM\tcount\n")
      for name in sorted([y.name for y in txome.transcripts],key=lambda x: calc.count(x), reverse=True):
         of4.write(name+"\t"+str(calc.TPM(name))+"\t"+str(calc.FPKM(name))+"\t"+str(calc.count(name))+"\n")
      of4.close()
      #if args.threads > 1: 
      #   p.close()
      #   p.join()

def prepare_gen(args,txome,rand,chunksize,cnt):
   for i in range(0,int(1+float(cnt)/float(chunksize))):
      yield (args,txome,rand.randint(0,10000000),chunksize)

def get_short_reads(vals):
   (args,txome,seed,chunk) = vals
   #fast forward some ammount
   """Emit the short reads first"""
   txe = TranscriptomeEmitter(txome,TranscriptomeEmitter.Options(seed=seed))
   if args.weights:
      weights = {}
      if args.weights[-3:]=='.gz': inf = gzip.open(args.weights)
      else: inf = open(args.weights)
      for line in inf:
         f = line.rstrip().split("\t")
         weights[f[0]] = float(f[1])
      txs = {}
      for tx in txome.transcripts: txs[tx.name] = tx.length
      for name in weights:
         weights[name] *= txs[name]
      txe.set_weights_by_dict(weights)
   else:
      weights = {}
      txs = {}
      for tx in txome.transcripts: txs[tx.name] = tx.length
      txe.set_weights_by_dict(weights)
   reademit = ReadEmitter(txe)
   shortreads = []
   sp = args.short_read_insert_size
   reademit.cutter.set_custom(sp[0],sp[1],sp[2])
   if args.short_read_error_rate: 
      emfr = ErrorMakerFlatRate(rate=args.short_read_error_rate,rand=reademit.options.rand)
      reademit.add_error_maker(emfr)
   for i in range(0,chunk):
      e = reademit.emit(args.short_read_length)
      shortreads.append(e)
   return shortreads

def get_long_reads(vals):
   (args,txome,seed,chunk) = vals
   txe = TranscriptomeEmitter(txome,TranscriptomeEmitter.Options(seed=seed))
   if args.weights:
      weights = {}
      if args.weights[-3:]=='.gz': inf = gzip.open(args.weights)
      else: inf = open(args.weights)
      for line in inf:
         f = line.rstrip().split("\t")
         weights[f[0]] = float(f[1])
      txe.set_weights_by_dict(weights)
   reademit = ReadEmitter(txe)
   longreads = []
   if args.long_read_error_rate: 
      emfr = ErrorMakerFlatRate(rate=args.long_read_error_rate,rand=reademit.options.rand)
      reademit.add_error_maker(emfr)
   lp = args.long_read_insert_size
   reademit.cutter.set_custom(lp[0],lp[1],lp[2])
   for i in range(0,chunk):
      e = reademit.emit(1)
      longreads.append(e)
   return longreads
      
def do_inputs():
  parser = argparse.ArgumentParser(description="simulate transcriptome reads",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  group1 = parser.add_argument_group('Required',description='You must set these')
  group1.add_argument('-o','--output',required=True,help="output folder")
  group1.add_argument('-t','--transcriptome',required=True,help="reference gpd transcriptome source")
  group1.add_argument('-g','--genome',required=True,help="reference genome fasta")
  group1.add_argument('--short_read_count',type=int,help="Number of PE short reads to generate")
  group1.add_argument('--long_read_count',type=int,help="Number of long reads to generate")
  group4 = parser.add_argument_group('Recommended',description='You may want to set these')
  group4.add_argument('--weights',help="set weights to choose transcripts from tsv. If not specified, all transcriptome transcripts are weighted evenly")
  group2 = parser.add_argument_group('ShortReads')
  group2.add_argument('--short_read_length',type=int,default=150,help="how long short reads should try to be (maybe shorter for short fragment)")
  group2.add_argument('--short_read_insert_size',nargs=3,type=float,default=[150,290,290],metavar=('MIN','MU','SIGMA'),help="minimum size and gaussian parameters mu and sigma")
  mgroup2 = group2.add_mutually_exclusive_group()
  mgroup2.add_argument('--short_read_error_rate',type=float,default=0,help="error rate for short reads")
  mgroup2.add_argument('--short_read_error_profile',nargs=2,metavar=('PROFILE','SCALE'),help="asssign an error profile from a file")
  mgroup2.add_argument('--short_read_specific_error',help="set a specific type of error")
  group3 = parser.add_argument_group('LongReads')
  group3.add_argument('--long_read_insert_size',nargs=3,type=float,default=[1000,4000,500],metavar=('MIN','MU','SIGMA'),help="minimum size and gaussian parameters mu and sigma")
  #group3.add_argument('--long_read_error_rate',type=float,default=0,help="General error rate to set for long reads")
  mgroup3 = group3.add_mutually_exclusive_group()
  mgroup3.add_argument('--long_read_error_rate',type=float,default=0,help="error rate for long reads")
  mgroup3.add_argument('--long_read_error_profile',nargs=2,metavar=('PROFILE','SCALE'),help="asssign an error profile from a file")
  mgroup3.add_argument('--long_read_specific_error',help="set a specific type of error")
  parser.add_argument('--threads',type=int,default=cpu_count(),help="Number of threads to convert names")
  parser.add_argument('--seed',type=int,help="set a seed by integer for a determinisitc set")
  args = parser.parse_args()
  if not args.short_read_count and not args.long_read_count:
     parser.error('--short_read_count or --long_read_count (or both) must be specified')
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
