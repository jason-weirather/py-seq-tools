"""simple unix sort for various formats

   sam and bam sort require samtools but reorders by coordinate instead of
   reordering by header order

"""
import argparse, sys, os, re
from shutil import rmtree
#from multiprocessing import cpu_count
from tempfile import mkdtemp, gettempdir
from subprocess import Popen, PIPE
from seqtools.format.sam import SAMStream, sort_header
from multiprocessing import cpu_count

def main(args):
  #do our inputs
  # Temporary working directory step 3 of 3 - Cleanup
  #sys.stderr.write("working in: "+args.tempdir+"\n")
  if args.bam or args.sam:
    do_sam(args)
    return
  cmd = "sort -T "+args.tempdir+'/'
  if args.psl:
    if args.name:
      cmd = "sort -k10,10 -T "+args.tempdir+'/'
    else:
      cmd = "sort -k14,14 -k15,15n -k16,16n -k9,9 -T "+args.tempdir+'/'
  if args.bed:
    cmd = "sort -k1,1 -k2,2n -T "+args.tempdir+'/'
  if args.gpd:
    if args.name:
      cmd = "sort -k1,1 -k2,2 -T "+args.tempdir+'/'
    else:
      cmd = "sort -k3,3 -k5,5n -k6,6n -k4,4 -T "+args.tempdir+'/'
  # Setup inputs 
  if args.input == '-':
    args.input = sys.stdin
  else:
    args.input = open(args.input)
  # Setup outputs
  if args.output:
    args.output = open(args.output,'w')
  else:
    args.output = sys.stdout
  p = Popen(cmd.split(),stdout=args.output,stdin=PIPE)
  for line in args.input:
    p.stdin.write(line)
  #p.stdin.close()
  #p.wait()
  p.communicate()
  if not args.specific_tempdir:
    rmtree(args.tempdir)

# special case for the sam type
def do_sam(args):
  if args.input != '-':
    m = re.search('\.bam$',args.input)
    if not m:  
      sys.stderr.write("ERROR input expects bam unless piping to stdin.. then SAM with header\n")
      sys.exit()
  of = sys.stdout
  cmdout2 = 'samtools view -S -h - '
  if args.output:
     of = open(args.output,'w')
  cmdout = 'samtools sort - '
  if args.threads:  cmdout += ' -@ '+str(args.threads)
  inf = None
  if args.input == '-':
    inf = sys.stdin
  else:
    cmd = 'samtools view -h '+args.input
    p = Popen(cmd.split(),stdout=PIPE,bufsize=1)
    inf = p.stdout
  s = SAMStream(inf)
  header = sort_header(s.header.text.rstrip())
  cmd2 = 'samtools view -Sb -'
  if args.output:
     if args.output[-4:] != '.bam': 
       pout2 = Popen(cmdout2.split(),stdin=PIPE,stdout=of)
       pout = Popen(cmdout.split(),stdin=PIPE,stdout=pout2.stdin)
     else:
       pout = Popen(cmdout.split(),stdin=PIPE,stdout=of)
  else:
     pout2 = Popen(cmdout2.split(),stdin=PIPE,stdout=of)
     pout = Popen(cmdout.split(),stdin=PIPE,stdout=pout2.stdin)
  p2 = Popen(cmd2.split(),stdin=PIPE,stdout=pout.stdin)
  p2.stdin.write(header.rstrip()+"\n")
  for sam in s:
    p2.stdin.write(str(sam)+"\n")
  p2.communicate()
  pout.communicate()
  if args.input != '-':
    p.communicate()
  return

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
  # Setup command line inputs
  parser=argparse.ArgumentParser(description="Sort these files by chromosome alphabetical, then start then end coordinate")
  parser.add_argument('input',help="INPUT FILE or '-' for STDIN")
  parser.add_argument('-o','--output',help="OUTPUTFILE or STDOUT if not set")
  parser.add_argument('--name',action='store_true',help="Sort by query name rather than location.  For GenePred this will default to gene name then the transcript name.")
  parser.add_argument('--threads',type=int,default=cpu_count(),help="INT number of threads to run. Default is system cpu count")
  # Temporary working directory step 1 of 3 - Definition
  group2 = parser.add_mutually_exclusive_group()
  group2.add_argument('--gpd',action='store_true')
  group2.add_argument('--bed',action='store_true')
  group2.add_argument('--psl',action='store_true')
  group2.add_argument('--bam',action='store_true',help="bam if file or sam if something else.")
  group2.add_argument('--sam',action='store_true',help="bam if file or sam if something else.")
  group = parser.add_mutually_exclusive_group()
  group.add_argument('--tempdir',default=gettempdir(),help="The temporary directory is made and destroyed here.")
  group.add_argument('--specific_tempdir',help="This temporary directory will be used, but will remain after executing.")
  args = parser.parse_args()

  # Temporary working directory step 2 of 3 - Creation
  setup_tempdir(args)
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

