#!/usr/bin/env python
from seqtools.format.gtf import GTF
import sys, argparse, gzip, re



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
  included = []
  sparse = []
  for line in inf:
     if re.match('#',line): continue
     gtf = GTF(line)
     vals = gtf.dict
     if args.feature and vals['feature'] != args.feature: continue
     attrib = vals['attribute']
     del vals['attribute']
     for k in vals:
       if k not in included: included.append(k)
     for k in attrib:
       knew = k
       if args.attribute_prefix: knew = args.attribute_prefix + k
       if knew not in included: included.append(k)
       if knew in vals: 
         raise ValueError('cant do a repeated value')
       vals[knew] = attrib[k]
     sparse.append(vals)
  if len(sparse)==0: return # nothing there?
  of.write("\t".join(included)+"\n")
  for line in sparse:
    v = [None if x not in line else line[x] for x in included]
    for value in [str(x) for x in v if x]: 
      if re.search('\t',value): raise ValueError('cant do internal tabs')
    of.write("\t".join(['' if x is None else str(x) for x in v])+"\n")
  of.close()


def do_inputs():
  # Setup command line inputs
  parser=argparse.ArgumentParser(description="Convert a gtf into a tsv",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="Use - for STDIN")
  parser.add_argument('--feature',help="Specify a feature to use from GTF to make into table")
  parser.add_argument('--output','-o',help="Specifiy path to write output or don't set for STDOUT")
  parser.add_argument('--attribute_prefix',help="Prefix attributes with this")
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
