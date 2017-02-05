import uuid, sys, time, re
import seqtools.structure.transcript
from seqtools.range import GenomicRange
from subprocess import Popen, PIPE
from collections import namedtuple

class GPD(seqtools.structure.transcript.Transcript):
  """ This whole format is a subclass of the Transcript subclass

  :param gpd_line:
  :type gpd_line: string
  """
  def __init__(self,gpd_line,options=None):
    if not options: options = {'sequence':None,
                               'ref':None,
                               'payload':None}
    # Only store the line and ID at first.  
    self._line = gpd_line.rstrip()
    m = re.match('[^\t]+\t[^\t]+\t([^\t]+)\t[^\t]+\t([^\t]+)\t([^\t]+)',gpd_line)
    self.entries = GPD._line_to_entry(self._line)

    exs = [GenomicRange(self.entries['chrom'], 
                        self.entries['exonStarts'][i]+1,
                        self.entries['exonEnds'][i]) for i in range(0,self.entries['exonCount'])]
    super(GPD,self).__init__(exs,{
      'direction':self.entries['strand'],
      'name':self.entries['name'],
      'gene_name':self.entries['gene_name'],
      'sequence':options['sequence'],
      'ref':options['ref'],
      'payload':options['payload']
    })

  def __str__(self):
    return self.get_gpd_line()  

  def get_gpd_line(self):
    """output the original gpd line
    Overrides Structure.Transcript"""
    return self._line

  def get_line(self):
    return self._line

  @staticmethod
  def _line_to_entry(line):
    f = line.rstrip().split("\t")
    starts = [int(x) for x in f[9].rstrip(",").split(",")]
    finishes =  [int(x) for x in f[10].rstrip(",").split(",")]
    return {
      'gene_name':f[0],
      'name':f[1],
      'chrom':f[2],
      'strand':f[3],
      'txStart':int(f[4]),
      'txEnd':int(f[5]),
      'cdsStart':int(f[6]),
      'cdsEnd':int(f[7]),
      'exonCount':int(f[8]),
      'exonStarts':starts,
      'exonEnds':finishes}

class GPDStream:
  """Iterate over GPD entries"""
  def __init__(self,fh):
    self.fh = fh

  def read_entry(self):
    ln = self.fh.readline()
    if not ln: return False
    gpd = GPD(ln)
    return gpd

  def __iter__(self):
    return self

  def next(self):
    r = self.read_entry()
    if not r:
      raise StopIteration
    else:
      return r

class SortedOutputFile:
  """a stream to write to for outputing a sorted file

  :param filename: output file
  :param type: how to sort (default location)
  :param tempdir:
  :type filename: string
  :type type: string
  :type tempdir: string
  """
  def __init__(self,filename,type='location',tempdir=None):
    if type not in ['location','name']:
      sys.stderr.write("ERROR: must be type location or name\n")
      sys.exit()
    self._gz = False
    self._fh = open(filename,'w')
    self._sh = None
    if filename[-3:] == '.gz':
      self._gz = True
    self._pipes  = []
    scmd = "sort -k1,1 -k2,2"
    if type == 'location':
      scmd = "sort -k3,3 -k5,5n -k6,6n -k4,4"
    if tempdir: scmd += " -T "+tempdir.rstrip('/')+'/'
    if self._gz:
      cmd1 = "gzip"
      p1 = Popen(cmd1.split(),stdout=self._fh,stdin=PIPE,close_fds=True)
      p2 = Popen(scmd.split(),stdout=p1.stdin,stdin=PIPE,close_fds=True)
      self._sh = p2.stdin
      self._pipes = [p2,p1]
    else:
      p = Popen(scmd.split(),stdout=self._fh,stdin=PIPE)
      self._sh = p.stdin
      self._pipes = [p]
  def write(self,value):
    self._sh.write(value)
  def close(self):
    #self._sh.flush()
    #self._sh.close()
    for p in self._pipes:
      #p.stdin.flush()
      #p.stdin.close()
      p.communicate()
    #self._pipes[0].stdin.flush()
    #self._pipes[0].stdin.close()
    #self._pipes[1].stdin.flush()
    #self._pipes[1].stdin.close()
    self._fh.close()
