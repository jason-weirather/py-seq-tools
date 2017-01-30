"""Convert a sam/bam to a psl, gzip output implied from output filename

   .. warning:: Requires a reference fasta now but we could change this behavior in the future to extract the necessary target lengths from the header


"""
import sys, argparse, inspect, os, gzip

from seqtools.format.fasta import FASTAData
from seqtools.format.sam import BAMFile, SAMStream

def main(args):
  of = sys.stdout
  if args.output:
    if args.output[-3:] == '.gz':
      of = gzip.open(args.output,'w')
    else: of = open(args.output,'w')

  ref = None
  if args.reference:
    ref = FASTAData(open(args.reference,'rb').read())
  
  if args.input == '-':
    args.input = SAMStream(sys.stdin,reference=ref)
  else: args.input = BAMFile(args.input,reference=ref)
  for e in args.input:
    if e.is_aligned():
      of.write(e.get_PSL()+"\n")
  of.close()


def do_inputs():
  parser = argparse.ArgumentParser(description="",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="Use - for STDIN or specify a BAM file")
  parser.add_argument('-r','--reference',help="Reference fasta",required=True)
  parser.add_argument('-o','--output',help="specify the output file or don't set for STDOUT")
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

