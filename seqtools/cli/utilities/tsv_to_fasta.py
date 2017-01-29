"""Undo the fasta_to_tsv command and put it back in fasta format"""
import argparse, sys, os, re

def main(args):
  inf = sys.stdin
  of = sys.stdout
  if args.input != '-':
    inf = open(args.input)
  if args.output:
    of = open(args.output,'w')
  for tsv in inf:
    f = tsv.rstrip().split("\t")
    if len(f) != 2:
      sys.stderr.write("only expecting two fields per line\n")
      sys.exit()
    of.write('>'+f[0]+"\n"+f[1]+"\n")
  of.close()

def do_inputs():
  # Setup command line inputs
  parser=argparse.ArgumentParser(description="Convert a tsv into a fasta assuming it is header <tab> sequence format",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
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
