"""Convert a sam/bam to a psl, gzip output implied from output filename

   .. warning:: Requires a reference fasta now but we could change this behavior in the future to extract the necessary target lengths from the header


"""
import sys, argparse, inspect, os, gzip

#from seqtools.format.fasta import FASTAData
from seqtools.format.sam import SAMStream
from seqtools.format.sam.bam.files import BAMFile

def main(args):
  of = sys.stdout
  if args.output:
    if args.output[-3:] == '.gz':
      of = gzip.open(args.output,'w')
    else: of = open(args.output,'w')

  if args.input == '-':
    args.input = SAMStream(sys.stdin)
  else: args.input = BAMFile(args.input)
  for e in args.input:
    if e.is_aligned():
      of.write(e.get_target_transcript(args.min_intron_size).get_gpd_line()+"\n")
  of.close()


def do_inputs():
  parser = argparse.ArgumentParser(description="",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="Use - for STDIN or specify a BAM file")
  parser.add_argument('-o','--output',help="specify the output file or don't set for STDOUT")
  parser.add_argument('--min_intron_size',type=int,default=68,help="the minimum intron size")
  args = parser.parse_args()
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

