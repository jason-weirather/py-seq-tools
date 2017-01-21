#!/usr/bin/env python
import argparse, sys, re
from seqtools.format.FASTQ import FASTQStream


def main(args):
  inf = sys.stdin
  of = sys.stdout
  if args.input != '-':
    inf = open(args.input)
  if args.output:
    of = open(args.output,'w')
  stream = FASTQStream(inf)
  for fq in stream:
    for l in fq.lines:
      if re.match('[\t]',l):
        sys.stderr.write('ERROR: there is a tab in the header\n'+l+"\n")
        sys.exit()
    of.write("\t".join(fq.lines)+"\n")
  of.close()

def do_inputs():
  # Setup command line inputs
  parser=argparse.ArgumentParser(description="Simply convert a fastq file to a tsv",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="Use - for STDIN")
  parser.add_argument('--output','-o',help="Specifiy path to write output")
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
