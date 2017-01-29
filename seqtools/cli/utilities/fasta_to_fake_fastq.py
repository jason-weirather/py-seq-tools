"""Convert a FASTA file into a FASTQ file. You can designate what to include
   in the quality score by setting the --ascii paramater (default 'I')"""
import argparse, sys
from seqtools.format.fasta import FASTAStream


def main(args):
  inf = sys.stdin
  of = sys.stdout
  if args.input != '-':
    inf = open(args.input)
  if args.output:
    of = open(args.output,'w')
  stream = FASTAStream(inf)
  for fa in stream:
    seq = fa.seq.replace("\n",'')
    of.write('@'+fa.header+"\n"+seq+"\n+\n"+args.ascii*len(seq)+"\n")
  of.close()

def do_inputs():
  # Setup command line inputs
  parser=argparse.ArgumentParser(description="Simply convert a fasta file to a fastq file with dummy quality",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="Use - for STDIN")
  parser.add_argument('--output','-o',help="Specifiy path to write output")
  parser.add_argument('--ascii',default='I',help="quality base to use")
  args = parser.parse_args()
  if len(args.ascii) != 1:
    sys.stderr.write('ERROR expecting just one character for ascii')
    sys.exit()
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
