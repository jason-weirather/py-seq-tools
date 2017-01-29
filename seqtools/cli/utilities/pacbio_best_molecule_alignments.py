"""Take alignments of PacBio data highest quality, any quality (low quality), 
or subreads and output the best aligned molecules

The purpose of this is to produce an output with only one read per molecule,
and that molecule representing the best alignment of that read.  Now the 'best'
alignment of a read is based on the following prioritized criteria.  Once a molecule is selected from one level, it will not be selected from another level.

Inputs are alignment in BAM format.

1. Aligned - The HQ (high quality) uncorrected alignmented reads are considered the best.  They are already single molecule. These are considered better than the HQ corrected alignments, because the HQ data should be set high enough that the correction process introduces more errors than it corrects, which would be expected with SNPs being artificially changed in correction at some point.
2. Aligned - The HQ corrected alignment is the next best.  If its available it certainly should be very high quality data, just not quite as good as the HQ uncorrected.
3. Aligned - The AQ (any quality) corrected is the next best.  These often include the HQ data, but since those have already been added if they are avaiable, its not an issue.  Any additional corrected alignments will be added.
4. Aligned - The AQ uncorrected reads are next.  These will be very diverse in their error rates.
5. Aligned - The corrected subreads are next.  The best aligned subread from each molecule is selected, (chosen by most bases aligned in a single path).
6. Aligned - The uncorrected subreads are the last group used.  Again they are chosen based on the best aligned subread from all the subreads of a molecule.
7. Unaligned - HQ - Any reads not aligned from the above sets that remain are added
8. Unaligned - HQ corrected - adding any more possible unaligned reads
9. Unaligned - AQ corrected - Any reads not aligned are added
10. Unaligned - AQ - Any reads not aligned are added
11. Unaligned - corrected subreads are added. since all will have an aligned base count of zero, the selection of which subread to add is based on subread length rather than aligned bases
12. Unaligned - uncorrected subreads are added. since all will have an aligned base count of zero, the selection of which subread to add is based on subread length rather than aligned bases

So basically we prioritize CCS reads so good they shouldn't be corrected first, then CCS reads that are corrected, then CCS reads, then best aligned subreads, and if thigns were not aligned we take the CCS read first, and the longest subread of the unaligned molecules last.

Unless a '.bam' file is specified the output is SAM format with the header derived from any one of the entered alignments.

"""
import argparse, sys, os, gzip, re
from seqtools.format.pacbio import PacBioReadName
from seqtools.format.sam import SAM
from subprocess import Popen, PIPE

_nameprog = re.compile('^(\S+)')

def main(args):
  of = sys.stdout
  if args.output and args.output[-4:] == '.bam':
    cmd = 'samtools view -Sb - -o '+args.output
    pof = Popen(cmd.split(),stdin=PIPE)
    of = pof.stdin
  elif args.output:
    of = open(args.output,'w')
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
  negative_filter = get_best_set(negative_filter,'-F 4',of,args,True)
  """After traversing all the aligned reads we can do it all over again
     this time with the unaligned portion of reads"""

  """ Finally go through the reads that did NOT ALIGN to get anything left"""
  get_best_set(negative_filter,'-f 4',of,args,False)
  if args.output and args.output[-4:] == '.bam':
    pof.communicate()
  else:
    of.close()

def get_best_set(negative_filter,flag,of,args,aligned):
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
    negative_filter = _do_subread_set(flag,args.subreadsCorrected,of,negative_filter,aligned)

  if args.subreads:
    """If we have subreads reads use those but we need to pick the best alignment for each molecule.
       The first pass through is just to get information on which is the best alignment"""
    negative_filter = _do_subread_set(flag,args.subreads,of,negative_filter,aligned)
  return negative_filter

def _do_subread_set(flag,input_file,of,negative_filter,aligned):
    best = {}
    cmd = 'samtools view '+flag+' '+input_file
    sys.stderr.write(cmd+"\n")
    p = Popen(cmd.split(),stdout=PIPE)
    z = 0
    for line in p.stdout:
      z += 1
      if z%10000==0: sys.stderr.write(str(z) + " subread alignment paths scanned for alignment length\r")
      pbname = PacBioReadName(_nameprog.match(line).group(1))
      mol = pbname.get_molecule()
      if mol in negative_filter: continue
      name = pbname.name()
      sam = SAM(line)
      c = 0 # aligned base count if we are aligned, subread length if we are not aligned
      if aligned:
        c = sam.get_aligned_bases_count()
      else:
        c = sam.get_query_length()
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
      if z%10000==0: sys.stderr.write(str(z) + " subreads alignment paths scanned during selected for best\r")
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
  parser.add_argument('--AQCorrected',help="Any quality CCS reads corrected. This is higher priority to use than the uncorrected any quality.")
  parser.add_argument('--subreads',help="the subreads")
  parser.add_argument('--subreadsCorrected',help="Any quality subread corrected.")
  parser.add_argument('--output','-o',help="Specifiy path to write index")
  args = parser.parse_args()
  if not (args.HQ or args.HQCorrected or args.AQ or args.AQCorrected \
       or args.subreads or args.subreadsCorrected):
    parser.error('must one or more alignment file --HQ --HQCorrected --AQ --AQCorrected --subreads --subreadsCorrected')
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
