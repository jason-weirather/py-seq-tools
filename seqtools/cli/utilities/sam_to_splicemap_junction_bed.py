#!/usr/bin/env python
"""Read any SAM/BAM and output a junction bed file in the splicemap format"""

import argparse, sys, re
from seqtools.format.sam import is_header, SAM
from seqtools.format.fasta import FASTAData
from subprocess import Popen, PIPE

def main(args):

  # get our reference genome
  sys.stderr.write("reading reference genome\n")

  g = FASTAData(open(args.reference_genome).read())
  sys.stderr.write("finished reading reference genome\n")

  inf = sys.stdin
  if args.infile != '-':
    if args.infile[-4:] == '.bam':
      cmd = 'samtools view -h '+args.infile
      pinf = Popen(cmd.split(),stdout=PIPE)
      inf = pinf.stdout
    else:
      inf = open(args.infile)

  read_mapping_count = {}
  junctions = {}
  sys.stderr.write("reading through sam file\n")
  zall = 0
  zn = 0
  for line in inf:
    if is_header(line): continue
    sam = SAM(line)
    if not sam.is_aligned(): continue

    chrom = sam.value('rname')
    if chrom =='*': continue
    if chrom not in g.keys():
      sys.stderr.write("WARNING: "+chrom+" not in reference, skipping\n")
      continue
    mate = 'U'
    if sam.check_flag(int('0x40',16)):
      mate = 'L'
    elif sam.check_flag(int('0x80',16)):
      mate = 'R'
    actual_read = sam.value('qname')+"\t"+mate
    if actual_read not in read_mapping_count:
      read_mapping_count[actual_read] = 0
    read_mapping_count[actual_read] += 1
    has_intron = 0
    start_loc = sam.value('pos')
    current_loc = start_loc
    bounds  = []
    cigar = [{'val':x[0],'op':x[1]} for x in sam.get_cigar()]
    for i in range(0,len(cigar)):
      # No action is necessary for H and S since they do not change the starting position on the reference
      ce = cigar[i]
      if ce['op'] == 'N' and ce['val'] >= args.min_intron_size:
        has_intron = 1
        lbound = current_loc # should be the intron start base index-1
        current_loc += ce['val']
        rbound = current_loc # should be the second exon start base index-1
        right_size = cigar[i+1]['val']
        bounds.append([lbound,rbound,right_size])
      elif ce['op'] == 'D':
        current_loc += ce['val']
      #elif re.match('[=XMSHP]',ce['op']):
      elif re.match('[=XM]',ce['op']):
        current_loc += ce['val']
      elif re.match('[P]',ce['op']):
        sys.stderr.write("ERROR: P padding not supported yet\n")
        sys.exit() 
    if has_intron == 0: continue # there are no splices to report here
    #print actual_read
    #print d['cigar']
    #print d
    #print start_loc
    #print bounds
    for bound in bounds:
      zall += 1
      intronflank = g[chrom][bound[0]-1:bound[0]+1].upper() + '-' + \
                    g[chrom][bound[1]-3:bound[1]-1].upper()
      strand = ''
      if is_canon(intronflank): # its a positive strand
        strand = '+'
      elif is_revcanon(intronflank): # its a negative strand
        strand = '-'
      else:
        # We can't deal with the non-canonical splice sorry
        zn += 1
        sys.stderr.write("WARNING skipping non-canonical splice ("+str(zn)+"/"+str(zall)+")\r")
        continue
      # If we are still in we have successfully found a splice
      out_chrom = chrom
      out_start = bound[0]-51
      out_end = bound[1]+49
      out_name = '*' # this will be done later
      out_score = 50
      out_strand = strand
      out_thickStart = out_start
      out_thickEnd = out_end
      out_rgb = '0,0,0'
      out_block_count = 2
      out_block_sizes = '50,50'
      out_block_starts = '0,'+str(bound[1]-bound[0]+50)
      bed = []
      bed.append(out_chrom)
      bed.append(str(out_start))
      bed.append(str(out_end))
      bed.append(out_name)
      bed.append(str(out_score))
      bed.append(out_strand)
      bed.append(str(out_thickStart))
      bed.append(str(out_thickEnd))
      bed.append(out_rgb)
      bed.append(str(out_block_count))
      bed.append(out_block_sizes)
      bed.append(out_block_starts)
      entry = "\t".join(bed)
      if entry not in junctions:
        junctions[entry] = {}
        junctions[entry]['reads'] = set()
        junctions[entry]['positions'] = set()
        junctions[entry]['right_sizes'] = set()
      junctions[entry]['reads'].add(actual_read)
      junctions[entry]['positions'].add(sam.value('pos'))
      junctions[entry]['right_sizes'].add(bound[2])
  sys.stderr.write("\n")
  sys.stderr.write("finished reading sam\n")
  of = sys.stdout
  if args.output:
    of = open(args.output,'w')
  if len(junctions) > 0: # if we have stuff lets print a header
    of.write("track\tname=junctions\tdescription=\"SpliceMap junctions\" itemRgb=\"On\"\n")
  for entry in junctions:
    nR = len(junctions[entry]['reads'])
    width = max(junctions[entry]['right_sizes'])-min(junctions[entry]['right_sizes'])
    nNR = len(junctions[entry]['positions'])
    nUR = 0
    nMR = 0
    for read in junctions[entry]['reads']:
      if read_mapping_count[read] == 1:
        nUR += 1
      elif read_mapping_count[read] > 1:
        nMR += 1
      else:
        sys.stderr.write("ERROR: nonsense read count\n")
        return
    name = '('+str(nR)+')['+str(width)+'_'+str(nNR)+']('+str(nUR)+'/'+str(nMR)+')'
    bed = entry.split("\t")
    bed[3] = name
    of.write("\t".join(bed)+"\n")    

  if args.infile != '-':
    if args.infile[-4:] == '.bam':
      pinf.communicate()
    else:
      inf.close()

def is_canon(input):
  v = set()
  v.add('GT-AG')
  v.add('GC-AG')
  v.add('AT-AC')
  if input in v: return True
  return False

def is_revcanon(input):
  v = set()
  v.add('CT-AC')
  v.add('CT-GC')
  v.add('GT-AT')
  if input in v: return True
  return False

def do_inputs():
  parser = argparse.ArgumentParser(description="Read a SAM/BAM file and output a bed file in the format of junction_color.bed",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('-o','--output',help='FILENAME is output')
  parser.add_argument('--min_intron_size',type=int,default=68,help='minimum intron size')
  parser.add_argument('infile',help='FILENAME SAM or BAM if .bam and "-" SAM for STDIN')
  parser.add_argument('-r','--reference_genome',required=True,help='REQUIRED: FILENAME of the reference genome')
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
