#!/usr/bin/env python
"""Convert a BAM/SAM file into one that is compatible with applications
that require a splicemap output SAM.  

**Requires samtools.**

This makes the IDP and IDP-fusion in versions that depend on splicemap
to accept more universal inputs.

"""
import sys, argparse
from seqtools.format.sam import SAMStream
from subprocess import PIPE, Popen

def main(args):
  inf = sys.stdin
  if args.input != '-':
    cmd = 'samtools view -h '+args.input
    pinf = Popen(cmd.split(),stdout=PIPE)
    inf = pinf.stdout
  of = sys.stdout

  if args.output:
    if args.output[-4:] == '.bam':
      cmd = 'samtools view -Sb - -o '+args.output
      pof = Popen(cmd.split(),stdin=PIPE)
      of = pof.stdin
    else: of = open(args.output,'w')

  stream = SAMStream(inf)
  of.write(stream.header_text)
  z = 0
  for sam in stream:
    z+=1
    if not sam.is_aligned(): continue # just keep aligned 
    if sam.entries.pos == 0: continue
    cigar = sam.cigar_array
    leftoffset = 0
    leftcigaroffset = 0
    if cigar[0][1] == 'H':
      if cigar[1][1] == 'S':
        leftoffset = cigar[1][0]
        leftcigaroffset = 2
    elif cigar[0][1] == 'S':
      leftoffset = cigar[0][0]
      leftcigaroffset = 1
    rightoffset = 0
    rightcigaroffset = 0
    if cigar[-1][1] == 'H':
      if cigar[-2][1] == 'S':
        rightoffset = cigar[-2][0]
        rightcigaroffset = 2
    elif cigar[-1][1] == 'S':
      rightoffset = cigar[-1][0]
      rightcigaroffset = 1
    origseq = sam.entries.seq
    seq = '*'
    if origseq != '*':
      seq = origseq[leftoffset:len(origseq)-rightoffset]
    qual = '*'
    if sam.entries.qual != '*':
      qual = sam.entries.qual[leftoffset:len(sam.entries.qual)-rightoffset]
    neocigar = cigar[leftcigaroffset:len(cigar)-rightcigaroffset]
    neocigarstring = ''.join([str(x[0])+x[1] for x in neocigar]).replace('N','D')
    flag = "0"
    if sam.check_flag(16):
      flag = "16"
    ostr = ''
    ostr += sam.entries.qname + "\t"
    ostr += flag + "\t"
    ostr += sam.entries.rname + "\t"
    ostr += str(sam.entries.pos) + "\t"
    ostr +=  "0" + "\t"
    ostr += neocigarstring + "\t"
    ostr += '*' + "\t"
    ostr += "0" + "\t"
    ostr += "0\t" if not sam.entries.tlen else str(sam.entries.tlen) + "\t"
    ostr += seq + "\t"
    ostr += qual
    of.write(ostr+"\n")
    if z%1000==0: sys.stderr.write("processed "+str(z)+" lines. At: "+sam.entries.rname + ':'+str(sam.entries.pos)+"           \r")

  sys.stderr.write("\n")
  if args.input != '-':
    pinf.communicate()
  if args.output:
    if args.output[-4:] == '.bam': pof.communicate()
    else: of.close()

def do_inputs():
  # Setup command line inputs
  parser=argparse.ArgumentParser(description="Generate our .bgi index for a bam file",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="Input BAM or - for STDIN SAM")
  parser.add_argument('--output','-o',help="Specifiy path to output, will be bam format if '.bam' is specified, STDOUT sam if not set")
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
