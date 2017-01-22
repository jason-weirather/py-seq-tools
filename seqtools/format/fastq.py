import re, sys
import seqtools.sequence

#Iterable Stream
class FASTQStream:
  def __init__(self,fh):
    self.fh = fh
  def __iter__(self):
    return self
  def next(self):
    v = self.get_entry()
    if not v:
      raise StopIteration
    else:
      return v
  def get_entry(self):
    line1 = self.fh.readline().rstrip()
    if not line1: return None
    line2 = self.fh.readline().rstrip()
    if not line2: return None
    line3 = self.fh.readline().rstrip()
    if not line3: return None
    line4 = self.fh.readline().rstrip()
    if not line4: return None
    return FASTQ("\n".join([line1,line2,line3,line4]))

class FASTQ(seqtools.sequence.Seq):
  # v is the lines with the first '@' removed
  def __init__(self,v):
    self.lines = v.rstrip().split("\n")
    self.header = self.lines[0][1:]
    self.name = re.match('(\S+)',self.header).group(1)
    self.seq = self.lines[1]
    self.qual = self.lines[3]
  def __getitem__(self,key):
    if isinstance(key,slice):
      newseq = self.seq[key.start:min(key.stop,len(self.seq))]
      newqual = self.qual[key.start:min(key.stop,len(self.seq))]
      return FASTQ('@'+self.header+"\n"+newseq+"\n"+self.lines[2]+"\n"+newqual)
    return {'name':self.name,'seq':self.seq,'qual':self.lines[3]}[key]
  def rc(self):
    return FASTQ('@'+self.header+"\n"+seqtools.sequence.rc(self.seq)+"\n"+self.lines[2]+"\n"+self.qual[::-1]+"\n")
  def copy(self):
    return FASTQ(self.FASTQ)
  def FASTQ(self):
    return "\n".join(self.lines)+"\n"
  def __str__(self):
    return self.FASTQ().rstrip()
