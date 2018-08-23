#!/usr/bin/env python
"""gpd_to_gtf.py 

   INPUT GENEPRED FORMAT:

   genePred format that associates a gene name with gene prediction information
   No header on input (may support in the future)
   No comments on input (may support in future)

   (
    string  geneName;           "Name of gene as it appears in Genome Browser."
    string  name;               "Name of gene"
    string  chrom;              "Chromosome name"
    char[1] strand;             "+ or - for strand"
    uint    txStart;            "Transcription start position"
    uint    txEnd;              "Transcription end position"
    uint    cdsStart;           "Coding region start"
    uint    cdsEnd;             "Coding region end"
    uint    exonCount;          "Number of exons"
    uint[exonCount] exonStarts; "Exon start positions"
    uint[exonCount] exonEnds;   "Exon end positions"
   )

   OUTPUT GTF FORMAT:

   genePred format that associates a gene name with gene prediction information
   No header on output (may support in the future)
   No comments on output (may support in future)
   Frame is not implemented on output (may support in future)
   Source is an optional 2nd input

   (
    string  chrom;           "Chromsome name."
    string  source;          "Source e.g. hg19"
    string  feature;         "exon"
    uint    start;           "Feature start position"
    uint    end;             "Feature end position"
    float   score;           ". (empty)"
    char[1] strand;          "+ or - for strand"
    char[1] frame;           ". (empty) for exons"
    string  attributes;      "gene/transcript/exon ids"
   )'


   Jason Weirather 20140317
   This script will take genePred files exons and convert them to gtf files
"""
from __future__ import print_function
import sys, gzip, argparse

source = '.'
def main(args):
   linenum = 0
   inf = sys.stdin
   if args.input != '-':
      if args.input[-3:] == '.gz': inf = gzip.open(args.input)
      else: inf = open(args.input) 
   of = sys.stdout
   if args.output:
     if args.output[-3:] == '.gz': inf = gzip.open(args.output,'w')
     else: of = open(args.output,'w')
   for line in inf:
      print(line)
      linenum+=1
      if line.startswith("#"): continue
      vals  = line.rstrip("\r\n").split("\t")
      gene = vals[0]
      txn = vals[1]
      chrom = vals[2]
      strand = vals[3]
      txStart = vals[4]
      txEnd = vals[5]
      cdsStart = vals[6]
      cdsEnd = vals[7]
      exonCount = vals[8]
      exonStarts = vals[9]
      exonEnds = vals[10]
      exonStarts = exonStarts.rstrip(",").split(",")
      exonEnds = exonEnds.rstrip(",").split(",")
      for i in range (0,len(exonStarts)):
         exonnumber = i+1
         exonid = txn+'.'+str(exonnumber)
         of.write(chrom + "\t" + source + "\texon\t" + str(int(exonStarts[i])+1) +"\t" + exonEnds[i] + "\t.\t" + strand + "\t.\t" + 'gene_id "'+gene+'"; transcript_id "'+txn+'"; exon_number "'+str(exonnumber)+'"; exon_id "'+exonid+'"; gene_name "'+gene+'";'+"\n")
   of.close()
   inf.close()

def do_inputs():
  # Setup command line inputs
  parser=argparse.ArgumentParser(description="Convert a genepred into a gtf",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="Use - for STDIN")
  parser.add_argument('--output','-o',help="Specifiy path to write output or don't set for STDOUT")
  parser.add_argument('--source',default='.',help="fill in source")
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

