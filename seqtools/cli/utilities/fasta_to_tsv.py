#!/usr/bin/env python
import argparse, sys, os, re
from seqtools.Format.Fasta import FastaHandle


def main(args):
  inf = sys.stdin
  of = sys.stdout
  if args.input != '-':
    inf = open(args.input)
  if args.output:
    of = open(args.output,'w')
  stream = FastaHandle(inf)
  for fa in stream:
    if re.match('[\t]',fa.header):
      sys.stderr.write("ERROR: tab in header cannot convert to tsv")
      sys.stderr.write("\n")
    of.write(fa.header+"\t"+fa.seq.replace("\n",'')+"\n")
  of.close()

def do_inputs():
  # Setup command line inputs
  parser=argparse.ArgumentParser(description="Simply convert a fasta into tsv output",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="Use - for STDIN")
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
