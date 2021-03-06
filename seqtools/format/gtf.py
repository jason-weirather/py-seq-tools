# GTFBasics.py
#
# A module for holding functions for handling psl files
# 
# These include:
# line_to_entry() - read a line from a psl file into a dictionary
#           also includes a new array of adjusted coordinates 
#           in the same coordiante system as a positive strand
#           for easier comparison of query alignments
#

import re, sys
from collections import OrderedDict

class GTFFile:
  def __init__(self,filehandle):
    self.genes = {}
    self.transcripts = {}
    for line in filehandle:
        if re.match('^#',line): continue
        f = line.rstrip().split("\t")
        if f[2] != 'exon': continue
        e = line_to_entry(line)
        if not e: continue
        if 'gene_id' not in e['attributes']:
          sys.stderr.write("WARNING no gene_id attribute found for "+line+"\n"+str(e)+"\n")
          continue
        if 'transcript_id' not in e['attributes']:
          sys.stderr.write("WARNING no gene_id attribute found for "+line+"\n")
          continue
        gene_id = e['attributes']['gene_id']
        transcript_id = e['attributes']['transcript_id']
        if gene_id not in self.genes:
          self.genes[gene_id] = []
        self.genes[gene_id].append(e)
        if transcript_id not in self.transcripts:
          self.transcripts[transcript_id] = []
        self.transcripts[transcript_id].append(e)
    return
  def write_genepred(self,filehandle):
    for transcript in self.transcripts:
      #filehandle.write(str(transcript)+"\n")
      tlist = self.transcripts[transcript]
      elist = {}
      gene = '.'
      strand = '.'
      chrom = '.'
      for t in tlist:
        if not t['gff'][2].lower() == 'exon':
          continue
        start = int(t['gff'][3])-1
        end = int(t['gff'][4])
        if start > end: 
          sys.stderr.write("ERROR start bigger than end\n")
          sys.exit()
        elist[start] = end
        gene = t['attributes']['gene_id']
        strand = t['gff'][6]
        chrom = t['gff'][0]
      #print strand
      #print transcript
      starts = sorted(elist,key=lambda k: elist[k])
      first = starts[0]
      last = elist[starts[len(starts)-1]]
      ostring = gene + "\t" + transcript + "\t" + chrom + "\t" + strand + "\t" + str(first) + "\t" \
              + str(last) + "\t" + str(first) + "\t" + str(last) + "\t" \
              + str(len(starts)) + "\t" \
              + ",".join([str(x) for x in starts]) + ",\t" \
              + ",".join([str(elist[x]) for x in starts])+","
      filehandle.write(ostring+"\n")

class GTF:
  def __init__(self,line=None):
    self.entry = None
    if line: self.entry = line_to_entry(line)
  @property
  def dict(self):
    d = OrderedDict()
    gff = self.entry['gff']
    d['seqname'] = None if gff[0] == '.' else gff[0]
    d['source'] = None if gff[1] == '.' else gff[1]
    d['feature'] = None if gff[2] == '.' else gff[2]
    d['start'] = None if gff[3] == '.' else int(gff[3])
    d['end'] = None if gff[4] == '.' else int(gff[4])
    d['score'] = None if gff[5] == '.' else gff[5]
    d['strand'] = None if gff[6] == '.' else gff[6]
    d['frame'] = None if gff[7] == '.' else int(gff[7])
    d['attribute'] = self.entry['attributes']
    return d

def line_to_entry(line):
  f = line.rstrip().split("\t")
  gff_fields = f[:8]
  preattributes = re.split('\s*;\s*',f[8])
  attributes = OrderedDict()
  for attribute in preattributes:
    m = re.search('(\S+)\s*["\']([^\'"]+)["\']',attribute)
    if m:  
      attributes[m.group(1)] = m.group(2)
  if len(attributes.keys()) > 0:
    entry = {}
    entry['gff'] = gff_fields
    entry['attributes'] = attributes
    return entry
  return None
