"""Filter genePred files on different attributes"""
import argparse, sys, re, gzip
from seqtools.format.gpd import GPD, GPDStream

def main(args):
  name_list = set()
  gene_name_list = set()
  if args.names:
    with open(args.names) as inf:
      for line in inf:
        f = line.rstrip().split("\t")
        name_list.add(f[0])
  if args.gene_names:
    with open(args.gene_names) as inf:
      for line in inf:
        f = line.rstrip().split("\t")
        gene_name_list.add(f[0])
  if args.junctions:
    juncs = set()
    fh = None
    if args.junctions[-3:] == '.gz':
      fh = gzip.open(args.junctions)
    else: fh = open(args.junctions)
    stream = GPDStream(fh)
    for gpd in stream:
      juncs.add(gpd.get_junction_string())
  inf = sys.stdin
  if args.input != '-':
    inf = open(args.input)
  for line in inf:
    if re.match('^#',line): continue
    is_good = True
    g = GPD(line.rstrip())
    tot = g.get_length()
    if args.min_length:
      if tot < args.min_length:
        is_good = False
    if args.max_length:
      if tot > args.max_length:
        is_good = False
    if args.names:
      if g.value('name') not in name_list:
        is_good = False
    if args.gene_names:
      if g.value('gene_name') not in args.gene_name_list:
        is_good = False
    if args.junctions:
      if g.get_junction_string() not in juncs:
        is_good = False
    # If we are still here we can print
    if not args.invert:
      if is_good: print line.rstrip()
    else:
      if not is_good: print line.rstrip()
  
def do_inputs():
  parser = argparse.ArgumentParser(description="Filter a genepred by transcript length")
  parser.add_argument('input',help="Input '-' for STDOUT")
  parser.add_argument('--min_length',type=int,help="Minimum transcript length")
  parser.add_argument('--max_length',type=int,help="Maximum transcript length")
  parser.add_argument('--names',help="filter on a name list")
  parser.add_argument('--gene_names',help="filter on a gene name list")
  parser.add_argument('--junctions',help="GPD file containing transcript with junctions to filter on")
  parser.add_argument('-v','--invert',action='store_true',help='Invert search result')
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
