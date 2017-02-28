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
      if args.threads > 1: p = Pool(processes=args.threads)
      of1 = gzip.open(args.output+'/SR_output_1.fq.gz','w')
      of2 = gzip.open(args.output+'/SR_output_2.fq.gz','w')
      of3 = gzip.open(args.output+'/SR_original_1.fa.gz','w')
      of4 = gzip.open(args.output+'/SR_original_2.fa.gz','w')
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
      sys.stderr.write("\n")
      of1.close()
      of2.close()
      of3.close()
      of4.close()
      sys.stderr.write("\n")
      #if args.threads > 1: 
      #   p.close()
      #   p.join()

   if args.long_read_count:
      gen = prepare_gen(args,txome,rand,1000,args.long_read_count)
      if args.threads > 1: p = Pool(processes=args.threads)
      of1 = gzip.open(args.output+'/LR_output.fq.gz','w')
      of2 = gzip.open(args.output+'/LR_original.fa.gz','w')
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
      sys.stderr.write("\n")
      of1.close()
      of2.close()
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
   reademit = ReadEmitter(txe)
   shortreads = []
   reademit.cutter.set_sr_cuts()
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
   reademit = ReadEmitter(txe)
   longreads = []
   if args.long_read_error_rate: 
      emfr = ErrorMakerFlatRate(rate=args.long_read_error_rate,rand=reademit.options.rand)
      reademit.add_error_maker(emfr)
   reademit.cutter.set_lr_cuts()
   for i in range(0,chunk):
      e = reademit.emit(1)
      longreads.append(e)
   return longreads
      
def do_inputs():
  parser = argparse.ArgumentParser(description="simulate transcriptome reads",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('-o','--output',required=True,help="output folder")
  parser.add_argument('--threads',type=int,default=cpu_count(),help="Number of threads to convert names")
  parser.add_argument('-t','--transcriptome',required=True,help="reference gpd transcriptome source")
  parser.add_argument('-g','--genome',required=True,help="reference genome fasta")
  parser.add_argument('--seed',type=int,help="set a seed by integer for a determinisitc set")
  parser.add_argument('--short_read_count',type=int,help="Number of PE short reads to generate")
  parser.add_argument('--short_read_length',type=int,default=150,help="how long short reads should try to be (maybe shorter for short fragment)")
  parser.add_argument('--short_read_error_rate',type=float,default=0,help="error rate for short reads")
  parser.add_argument('--long_read_count',type=int,help="Number of long reads to generate")
  parser.add_argument('--long_read_error_rate',type=float,default=0,help="General error rate to set for long reads")
  args = parser.parse_args()
  return args
  
def external_cmd(cmd):
  cache_argv = sys.argv
  sys.argv = cmd.split()
  args = do_inputs()
  main(args)
  sys.argv = cache_argv

if __name__=="__main__":
  args = do_inputs()
  main(args)
