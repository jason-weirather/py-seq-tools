#!/usr/bin/env python
import argparse, sys, os, re
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
    if not args.inv:
      if args.left:
        fa = fa[args.left:]
      if args.right:
        fa = fa[:-1*args.right]
    else:
      if args.left:
        fa = fa[0:args.left]
      if args.right:
        fa = fa[-1*args.right:]
    of.write(fa.FASTA())
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
