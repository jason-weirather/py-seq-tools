#!/usr/bin/env python
import sys, gzip, argparse

#pre: a genepred file, a reference fasta file, an output fasta, optionally 'directionless' will reverse compliment entries on the negative strand
#post: writes to the output fasta the sequence of the transcripts in the genepred
#      produces transcripts in the proper orientation according to their direction in the genepred
#modifies fileIO

from seqtools.format.gpd import GPD
from seqtools.format.fasta import FASTAData

def main(args):
   inf = sys.stdin
   if args.input != '-':
      if args.input[-3:] == '.gz': inf = gzip.open(args.input)
      else: inf = open(args.input)
   of = sys.stdout
   if args.output:
      if args.output[-3:] =='.gz': of = gzip.open(args.output,'w')
      else: of = open(args.output,'w')
   fasta = None
   if args.genome[-3:] == '.gz': fasta = FASTAData(gzip.open(args.genome).read())
   else: fasta = FASTAData(open(args.genome).read())
   for line in inf:
      gpd = GPD(line)
      gpd.set_reference(fasta)
      ostr = '>'+gpd.transcript_name+"\n"
      ostr += str(gpd.sequence)+"\n"
      of.write(ostr)

def do_inputs():
   parser = argparse.ArgumentParser(description="Filter a genepred by transcript length")
   parser.add_argument('input',help="Input '-' for STDOUT")
   parser.add_argument('-g','--genome',required=True,help='Reference Genome Fasta, can be gzipped')
   parser.add_argument('-o','--output',help='Output file or STDOUT if not set')
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
