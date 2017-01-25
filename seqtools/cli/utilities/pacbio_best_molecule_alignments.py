"""Take alignments of PacBio data highest quality, any quality (low quality), 
   or subreads and output the best aligned molecules"""
import argparse, sys, os, gzip, re
from shutil import rmtree, copy
from multiprocessing import cpu_count, Pool
from tempfile import mkdtemp, gettempdir
from seqtools.format.pacbio import PacBioReadName
from seqtools.format.sam import SAM
from subprocess import Popen, PIPE

_nameprog = re.compile('^(\S+)')

def main(args):
  of = sys.stdout
  if args.output[-4:] == '.bam':
    cmd = 'samtools view -Sb - -o '+args.output
    pof = Popen(cmd.split(),stdin=PIPE)
    of = pof.stdin
  """Use the valid input file to get the header information."""

  header = None
  if args.HQ:
    cmd = 'samtools view -H '+args.HQ
    sys.stderr.write(cmd+"\n")
    header = Popen(cmd.split(),stdout=PIPE).communicate()[0]
    of.write(header)
  if (not header) and args.HQCorrected:
    cmd = 'samtools view -H '+args.HQCorrected
    sys.stderr.write(cmd+"\n")
    header = Popen(cmd.split(),stdout=PIPE).communicate()[0]
    of.write(header)
  if (not header) and args.AQ:
    cmd = 'samtools view -H '+args.AQ
    sys.stderr.write(cmd+"\n")
    header = Popen(cmd.split(),stdout=PIPE).communicate()[0]
    of.write(header)
  if (not header) and args.AQCorrected:
    cmd = 'samtools view -H '+args.AQCorrected
    sys.stderr.write(cmd+"\n")
    header = Popen(cmd.split(),stdout=PIPE).communicate()[0]
    of.write(header)
  if (not header) and args.subreads:
    cmd = 'samtools view -H '+args.subreads
    sys.stderr.write(cmd+"\n")
    header = Popen(cmd.split(),stdout=PIPE).communicate()[0]
    of.write(header)
  if (not header) and args.subreadsCorrected:
    cmd = 'samtools view -H '+args.subreadsCorrected
    sys.stderr.write(cmd+"\n")
    header = Popen(cmd.split(),stdout=PIPE).communicate()[0]
    of.write(header)

  _nameprog = re.compile('^(\S+)')
  negative_filter = set() # remove these
  """ Next read throught he alignments THAT ALIGNED in order of priority"""
  negative_filter = get_best_set(negative_filter,'-F 4',of,args)

  """ Finally go through the reads that did NOT ALIGN to get anything left"""
  get_best_set(negative_filter,'-f 4',of,args)
  if args.output[-4:] == '.bam':
    pof.communicate()
  else:
    of.close()

def get_best_set(negative_filter,flag,of,args):
  if args.HQ:
    """If we have the highest quality ccs reads use those first"""
    cmd = 'samtools view '+flag+' '+args.HQ
    sys.stderr.write(cmd+"\n")
    p = Popen(cmd.split(),stdout=PIPE)
    negative_filter = _traverse_unobserved(p.stdout,negative_filter,of)
    p.communicate()
    sys.stderr.write("molecules written: "+str(len(negative_filter))+"\n")
  if args.HQCorrected:
    """If we have the highest quality corrected ccs reads use those next"""
    cmd = 'samtools view '+flag+' '+args.HQCorrected
    sys.stderr.write(cmd+"\n")
    p = Popen(cmd.split(),stdout=PIPE)
    negative_filter = _traverse_unobserved(p.stdout,negative_filter,of)
    p.communicate()
    sys.stderr.write("molecules written: "+str(len(negative_filter))+"\n")
  if args.AQCorrected:
    """If we have the any quality corrected ccs reads use those next"""
    cmd = 'samtools view '+flag+' '+args.AQCorrected
    sys.stderr.write(cmd+"\n")
    p = Popen(cmd.split(),stdout=PIPE)
    negative_filter = _traverse_unobserved(p.stdout,negative_filter,of)
    p.communicate()
    sys.stderr.write("molecules written: "+str(len(negative_filter))+"\n")
  if args.AQ:
    """If we have the lower quality ccs reads use those next"""
    cmd = 'samtools view '+flag+' '+args.AQ
    sys.stderr.write(cmd+"\n")
    p = Popen(cmd.split(),stdout=PIPE)
    negative_filter = _traverse_unobserved(p.stdout,negative_filter,of)
    p.communicate()
    sys.stderr.write("molecules written: "+str(len(negative_filter))+"\n")

  if args.subreadsCorrected:
    """If we have corrected subreads reads use those but we need to pick the best alignment for each molecule.
       The first pass through is just to get information on which is the best alignment"""
    negative_filter = _do_subread_set(flag,args.subreadsCorrected,of,negative_filter)

  if args.subreads:
    """If we have subreads reads use those but we need to pick the best alignment for each molecule.
       The first pass through is just to get information on which is the best alignment"""
    negative_filter = _do_subread_set(flag,args.subreads,of,negative_filter)
  return negative_filter

def _do_subread_set(flag,input_file,of,negative_filter)
    best = {}
    cmd = 'samtools view '+flag+' '+input_file
    sys.stderr.write(cmd+"\n")
    p = Popen(cmd.split(),stdout=PIPE)
    z = 0
    for line in p.stdout:
      z += 1
      if z%1000==0: sys.stderr.write(str(z) + " subread alignment paths scanned for alignment length\r")
      pbname = PacBioReadName(_nameprog.match(line).group(1))
      mol = pbname.get_molecule()
      name = pbname.name()
      if mol in negative_filter: continue
      sam = SAM(line)
      c = 0 # aligned base cout
      if sam.is_aligned():
        c = sam.get_aligned_bases_count()
      if mol not in best:
        best[mol] = [name,c]
      if c > best[mol][1]: best[mol] = [name,c]
    p.communicate()
    sys.stderr.write("\n")
    sys.stderr.write("Finished analyzing subread lengths\nWriting aligned subreads\n")
    """After getting all the best alignment counts we can traverse again
       to keep the best"""
    cmd = 'samtools view '+flag+' '+input_file
    sys.stderr.write(cmd+"\n")
    z = 0
    p = Popen(cmd.split(),stdout=PIPE)
    for line in p.stdout:
      z += 1
      if z%1000==0: sys.stderr.write(str(z) + " subreads alignment paths scanned during selected for best\r")
      pbname = PacBioReadName(_nameprog.match(line).group(1))
      mol = pbname.get_molecule()
      name = pbname.name()
      if mol in negative_filter: continue
      if not best[mol][0] == name: continue
      of.write(line)
    p.communicate()
    for mol in best: negative_filter.add(mol)
    sys.stderr.write("\n")
    sys.stderr.write("molecules written: "+str(len(negative_filter))+"\n")
  """After traversing all the aligned reads we can do it all over again
     this time with the unaligned portion of reads"""
  return negative_filter

def _traverse_unobserved(stream,negative_filter,of):
  """Go through a stream and print out anything not in observed set"""
  observed = set()
  for line in stream:
    name = PacBioReadName(_nameprog.match(line).group(1))
    if name.get_molecule() not in negative_filter: of.write(line)
    observed.add(name.get_molecule())
  return negative_filter|observed

def do_inputs():
  # Setup command line inputs
  parser=argparse.ArgumentParser(description="Get a set of best molecule alignments. Prioritize on HQ.  This script does not account for corrected reads.  That will be another script. Requires samtools.",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('--HQ',help="highest quality CCS reads. This is the highest priority to use")
  parser.add_argument('--HQCorrected',help="highest quality corrected CCS reads.  This is lower priority to use than the uncorrected high quality.")
  parser.add_argument('--AQ',help="any quality CCS reads of a broad range of qualities (usually would includes HQ). This is lower priority to use than the corrected any quality.")
  parser.add_argument('--subreads',help="the subreads")
  parser.add_argument('--AQCorrected',help="Any quality CCS reads corrected. This is higher priority to use than the uncorrected any quality.")
  parser.add_argument('--subreadsCorrected',help="Any quality subread corrected.")
  parser.add_argument('--output','-o',help="Specifiy path to write index")
  parser.add_argument('--threads',type=int,default=cpu_count(),help="INT number of threads to run. Default is system cpu count")
  # Temporary working directory step 1 of 3 - Definition
  group = parser.add_mutually_exclusive_group()
  group.add_argument('--tempdir',default=gettempdir(),help="The temporary directory is made and destroyed here.")
  group.add_argument('--specific_tempdir',help="This temporary directory will be used, but will remain after executing.")
  group.add_argument('--no_primary_search',action='store_true',help="Dont try to assign a primary flag")
  args = parser.parse_args()
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
