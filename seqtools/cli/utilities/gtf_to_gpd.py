#!/usr/bin/env python
import seqtools.format.gtf
import sys, argparse, gzip



# Pre:  A GTF filename, for a file with 'exon' features, and 'gene_id' and 'transcript_id' attributes.
# Post: Prints to stdout a genepred with transcripts

def main(args):
  of = sys.stdout
  inf = sys.stdin
  if args.input != '-': 
     if args.input[-3:] == '.gz':
       inf = gzip.open(args.input)
     else: 
       inf = open(args.input)
  if args.output: 
     if args.output[-3:] == '.gz':
       of = gzip.open(args.output,'w')
     else:
       of = open(args.output,'w')
  gtf = seqtools.format.gtf.GTFFile(inf)
  gtf.write_genepred(of)

def do_inputs():
  # Setup command line inputs
  parser=argparse.ArgumentParser(description="Convert a gtf into a gpd",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="Use - for STDIN")
  parser.add_argument('--output','-o',help="Specifiy path to write output or don't set for STDOUT")
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
