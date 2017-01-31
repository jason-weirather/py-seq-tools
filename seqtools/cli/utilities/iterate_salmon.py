"""Iteratively execute salmon until all transcripts with
   leave-one-out-ish until all transcripts meet expression criteria

   a nice addition (or option) would be to have the iteration 
   continuously testing the gene group to see if it should be 
   split apart more because of no linkage
"""
import argparse, sys, os, gzip
from shutil import rmtree, copy
from multiprocessing import cpu_count, Pool
from tempfile import mkdtemp, gettempdir
from seqtools.format.fasta import FASTAData
from seqtools.format.gpd import GPDStream
from seqtools.cli.utilities.fastq_to_salmon_quant import external_cmd as salmon
from subprocess import Popen, PIPE

def main(args):

  cmd_fastq = ' -p '+str(args.numThreads)
  if args.rU: cmd_fastq += ' --rU '+args.rU
  if args.r1: cmd_fastq += ' --r1 '+args.r1
  if args.r2: cmd_fastq += ' --r2 '+args.r2
  cmd_genome = ' '+args.genome
  z = 0
  current_gpd = args.gpd
  prev_size = -1
  while True:
    z += 1
    result = args.tempdir+'/result.'+str(z)
    # run an iteration
    cmd_gpd = ' '+current_gpd
    cmd = 'fastq_to_salmon_quant '+cmd_fastq+cmd_gpd+cmd_genome+' -o '+result
    sys.stderr.write(cmd)
    salmon(cmd)
    #check the result
    fail = set()
    expression = {}
    with open(result) as inf:
      header = inf.readline()
      for line in inf:
        f = line.rstrip().split("\t")
        tpm = float(f[4])
        gtpm = float(f[6])
        gene = f[0]
        transcript = f[1]
        expression[transcript] = tpm
        if tpm < args.minTPM:
          fail.add(transcript)
        if gtpm > 0 and tpm/gtpm < args.minFraction:
          fail.add(transcript)
    # now check the gpd
    genes = {}
    with open(current_gpd) as inf:
      for line in inf:
        f = line.rstrip().split("\t")    
        if f[0] not in genes: genes[f[0]] = []
        genes[f[0]].append(f)
    current_gpd = args.tempdir+'/ref.'+str(z)+'.gpd'
    tof = open(current_gpd,'w')
    # populated all the genes
    cnt = 0
    for gene in genes:
      expressed = [x for x in genes[gene] if x[1] in expression and expression[x[1]] > 0.0001]
      txs = sorted(expressed, key=lambda x: expression[x[1]])
      #print [expression[x[1]] for x in txs]
      failures = len([x for x in txs if x[1] in fail])
      dcount = max(1,int(float(failures)*args.decimate))
      for i in range(0,dcount):
        if len(txs) > 0:
          if txs[0][1] in fail: txs.pop(0)
      for tx in txs:
        tof.write("\t".join(tx)+"\n")
        cnt += 1
    tof.close()
    if cnt == 0 or cnt == prev_size: break
    prev_size = cnt
  # Now can output the final results
  of = sys.stdout
  if args.output: of = open(args.output,'w')
  with open(args.tempdir+'/result.'+str(z)) as inf:
    for line in inf: of.write(line)
  of.close()
  # Temporary working directory step 3 of 3 - Cleanup
  if not args.specific_tempdir:
    rmtree(args.tempdir)

def do_inputs():
  # Setup command line inputs
  parser=argparse.ArgumentParser(description="Quantify a genepred defined transcriptome using salmon",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('--rU',help="INPUT FASTQ/FASTA can be gzipped")
  parser.add_argument('--r1',help="INPUT Pair1 FASTQ/FASTA can be gzipped")
  parser.add_argument('--r2',help="INPUT Pair2 FASTQ/FASTA can be gzipped")
  parser.add_argument('gpd',help="transcriptome genepred")
  parser.add_argument('genome',help="reference fasta")
  parser.add_argument('--output','-o',required=True,help="Specifiy path write output or - for STDOUT")
  parser.add_argument('-p','--numThreads',type=int,default=cpu_count(),help="INT number of threads to run. Default is system cpu count")
  parser.add_argument('-t','--minTPM',type=float,default=1.0,help='smallest tpm for a transcript in final output')
  parser.add_argument('-f','--minFraction',type=float,default=0.1,help='smallest fraction of gene for a transcript to be in the end')
  parser.add_argument('--decimate',type=float,default=0.1,help='if number of failures is large decmiate it by this fraction')
  # Temporary working directory step 1 of 3 - Definition
  group = parser.add_mutually_exclusive_group()
  group.add_argument('--tempdir',default=gettempdir(),help="The temporary directory is made and destroyed here.")
  group.add_argument('--specific_tempdir',help="This temporary directory will be used, but will remain after executing.")
  args = parser.parse_args()
  if not ((args.r1 and args.r2) or args.rU):
    parser.error("Either rU or both r1 and r2 need to be specified")
  # Temporary working directory step 2 of 3 - Creation
  setup_tempdir(args)
  return args

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

def external_cmd(cmd):
  cache_argv = sys.argv
  sys.argv = cmd.split()
  args = do_inputs()
  main(args)
  sys.argv = cache_argv

if __name__=="__main__":
  args = do_inputs()
  main(args)
