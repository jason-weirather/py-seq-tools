"""Requires sorted gpd. Change names in format to a random assignment. Cluster locus and random assign a gene name"""
import argparse, sys, re, uuid, gzip
from seqtools.format.gpd import GPDStream
from seqtools.stream import LocusStream

def main(args):
  inf = sys.stdin
  of = sys.stdout
  if args.input != '-':
    if args.input[-3:] == '.gz': inf = gzip.open(args.input)
    else: inf = open(args.input)
  if args.output:
    of = open(args.output,'w')
  stream = LocusStream(GPDStream(inf))
  for rng in stream:
    gpds = rng.get_payload()
    overlap_groups = []
    while len(gpds) > 0:
      check = gpds.pop(0)
      # see if check overlaps any groups
      overlaps_one = False
      for i in range(0,len(overlap_groups)):
        group = overlap_groups[i]
        for member in group:
          if check.overlap_size(member) > 0:
            overlaps_one = True
            overlap_groups[i].append(check)
            break
        if overlaps_one: break
      if not overlaps_one:
        overlap_groups.append([check]) # add a group
    # now we have groups that overlap
    for group in overlap_groups:
      gname = str(uuid.uuid4())
      for member in group:
        f = str(member).rstrip().split("\t")
        f[0] = gname
        f[1] = str(uuid.uuid4())
        of.write("\t".join(f)+"\n")
    #f = line.rstrip().split("\t")
    #f[1] = str(uuid.uuid4())
    #of.write("\t".join(f)+"\n")
  of.close()

def do_inputs():
  # Setup command line inputs
  parser=argparse.ArgumentParser(description="Convert a tsv into a fastq",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
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
