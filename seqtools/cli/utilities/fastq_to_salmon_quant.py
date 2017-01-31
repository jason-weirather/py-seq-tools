"""Take a fastq/fasta and make a transcriptome quantification

   **Requires salmon**

   Pre: Requires

   1. A FASTQ/FASTA input file (can be gzipped)
   2. A transcriptome genepred
   3. A genome fasta

   See -h for detailed input description

   Post: Returns

   Writes an output table that has header and the following fields:

   1. geneName
   2. transcriptName
   3. length (transcript)
   4. EffectiveLength (transcript)
   5. TPM (transcript)
   6. NumReads(transcript)
   7. TPM (gene)
"""
import argparse, sys, os, gzip
from shutil import rmtree, copy
from multiprocessing import cpu_count, Pool
from tempfile import mkdtemp, gettempdir
from seqtools.format.fasta import FASTAData
from seqtools.format.gpd import GPDStream
from subprocess import Popen, PIPE

def main(args):

  """ First read the genome """
  sys.stderr.write("reading reference genome\n")
  ref = FASTAData(open(args.genome).read())
  sys.stderr.write("read in "+str(len(ref.keys()))+" chromosomes\n")

  """ Next make the transcriptome """
  txome = {}
  sys.stderr.write("write the transcriptome\n")
  inf = None
  if args.gpd[-3:] == '.gz':
    inf = gzip.open(args.gpd)
  else:
    inf = open(args.gpd)
  stream = GPDStream(inf)
  tof = open(args.tempdir+'/transcriptome.fa','w')
  z = 0
  for gpd in stream:
    z += 1
    if gpd.get_transcript_name() in txome:
      sys.stderr.write("WARNING already have a transcript "+gpd.get_transcript_name()+" ignoring line "+str(z)+" of the gpd\n")
      continue
    txome[gpd.get_transcript_name()] = gpd.get_gene_name()
    tof.write('>'+gpd.get_transcript_name()+"\n"+gpd.get_sequence(ref)+"\n")
  tof.close()
  inf.close()
  sys.stderr.write("wrote "+str(len(txome.keys()))+" transcripts\n")

  """Build the salmon index"""
  sys.stderr.write("building a salmon index\n")
  cmd = 'salmon index -p '+str(args.numThreads)+' -t '+args.tempdir+'/transcriptome.fa -i '+args.tempdir+'/salmon_index'
  p = Popen(cmd.split())
  p.communicate()
  sys.stderr.write("finished building the index\n")

  """Use the index to quanitfy"""
  sys.stderr.write("quanitfy reads\n")
  reads = ''
  if args.rU:
    reads = '-r '+args.rU
  else:
    reads = '-1 '+args.r1+' -2 '+args.r2
  cmd = 'salmon quant -p '+str(args.numThreads)+' -i '+args.tempdir+'/salmon_index -l A '+reads+' -o '+args.tempdir+'/output_quant'
  p = Popen(cmd.split())
  p.communicate()
  sys.stderr.write("finished quanitfying\n")

  """Now parse the salmon output to add gene name"""
  salmon = {}
  with open(args.tempdir+'/output_quant/quant.sf') as inf:
    header = inf.readline()
    for line in inf:
      f = line.rstrip().split("\t")
      # by each transcript name hold a data strcture of the information
      salmon[f[0]] = {'name':f[0],'length':int(f[1]),'EffectiveLength':float(f[2]),'TPM':float(f[3]),'NumReads':float(f[4])}
  genes = {}
  for name in salmon:
    gene = txome[name]
    if gene not in genes: genes[gene] = []
    genes[gene].append(salmon[name])
  genetot = {}
  for gene in genes:
    tot = sum([x['TPM'] for x in genes[gene]])
    genetot[gene] = tot
  ordered_gene_names = sorted(genetot.keys(), key=lambda x: genetot[x],reverse=True)
      
  """Collected enough information to make output"""
  sys.stderr.write("generating output\n")
  of = sys.stdout
  if args.output != '-':
    of = open(args.output,'w')
  of.write("geneName\ttranscriptName\tlength\tEffectiveLength\ttxTPM\tNumReads\tgeneTPM\n")
  for gene in ordered_gene_names:
    txs = sorted(genes[gene],key=lambda x: x['TPM'],reverse=True)
    for tx in txs:
      of.write(gene+"\t"+tx['name']+"\t"+str(tx['length'])+"\t"+str(tx['EffectiveLength'])+"\t"+str(tx['TPM'])+"\t"+str(tx['NumReads'])+"\t"+str(genetot[gene])+"\n")
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
  parser.add_argument('--output','-o',required=True,help="Specifiy path to write output or - for stdout")

  parser.add_argument('-p','--numThreads',type=int,default=cpu_count(),help="INT number of threads to run. Default is system cpu count")

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
