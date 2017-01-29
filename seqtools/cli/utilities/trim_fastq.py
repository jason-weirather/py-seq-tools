"""Trim a FASTQ file/stream ends of all entries (option to invert and only keep the ends)"""
import argparse, sys, os, re
from seqtools.format.fastq import FASTQStream

def main(args):
  inf = sys.stdin
  of = sys.stdout
  if args.input != '-':
    inf = open(args.input)
  if args.output:
    of = open(args.output,'w')
  stream = FASTQStream(inf)
  for fq in stream:
    if not args.inv:
      if args.left:
        fq = fq[args.left:]
      if args.right:
        fq = fq[:-1*args.right]
    else:
      if args.left:
        fq = fq[0:args.left]
      if args.right:
        fq = fq[-1*args.right:]
    of.write(fq.FASTQ())
  of.close()

def do_inputs():
  # Setup command line inputs
  parser=argparse.ArgumentParser(description="Trim bases off ends (or keep)",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="Use - for STDIN")
  group = parser.add_mutually_exclusive_group()
  group.add_argument('--left',type=int)
  group.add_argument('--right',type=int)
  parser.add_argument('--inv',action='store_true')
  parser.add_argument('--output','-o',help="Specifiy path to write index")
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
