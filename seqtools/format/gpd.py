import uuid, sys, time, re, os
import seqtools.structure.transcript
from seqtools.range import GenomicRange
from subprocess import Popen, PIPE
from collections import namedtuple


"""GPD options must be declared at the top level of the module to
   prevent problems with pickeling"""
GPDOptions = namedtuple('GPDOptions',
   ['sequence',
    'ref',
    'payload'])

"""GPDFields namedtuple must be defined at the top level to
   facilitate serialization"""
GPDFields = namedtuple('GPDFields',
   ['gene_name',
    'name',
    'chrom',
    'strand',
    'txStart',
    'txEnd',
    'cdsStart',
    'cdsEnd',
    'exonCount',
    'exonStarts',
    'exonEnds'])

class GPD(seqtools.structure.transcript.Transcript):
  """ This whole format is a subclass of the Transcript subclass

  :param gpd_line:
  :type gpd_line: string
  """
  def __init__(self,gpd_line,options=None):
    if not options: options = GPD.Options()
    self._options = options
    # Only store the line and ID at first.
    self._line = gpd_line.rstrip()
    m = re.match('[^\t]+\t[^\t]+\t([^\t]+)\t[^\t]+\t([^\t]+)\t([^\t]+)',gpd_line)
    self.entries = GPD._line_to_entry(self._line)

    exs = [GenomicRange(self.entries.chrom,
                        self.entries.exonStarts[i]+1,
                        self.entries.exonEnds[i]) for i in range(0,self.entries.exonCount)]
    super(GPD,self).__init__(exs,
                             super(GPD,self).Options(
      direction = self.entries.strand,
      name = self.entries.name,
      gene_name = self.entries.gene_name,
      sequence = options.sequence,
      ref = options.ref,
      payload = options.payload)
    )

  @staticmethod
  def Options(**kwargs):
     """ A method for declaring options for the class"""
     construct = GPDOptions #IMPORTANT!  Set this
     names = construct._fields
     d = {}
     for name in names: d[name] = None #default values
     for k,v in kwargs.items():
       if k in names: d[k] = v
       else: raise ValueError('Error '+k+' is not a property of these options')
     """Create a set of options based on the inputs"""
     return construct(**d)

  @property
  def name(self):
    return self.entries.name

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
    return GPDFields(
      f[0],
      f[1],
      f[2],
      f[3],
      int(f[4]),
      int(f[5]),
      int(f[6]),
      int(f[7]),
      int(f[8]),
      starts,
      finishes)

class GPDStream:
  """Iterate over GPD entries"""
  def __init__(self,fh):
    self.fh = fh

  def read_entry(self):
    ln = self.fh.readline()
    if not ln: return False
    gpd = GPD(ln.decode('utf-8'))
    return gpd

  def __iter__(self):
    return self

  def __next__(self):
    return self.next()

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
    scmd = ['sort','-k1,1','-k2,2']
    if type == 'location':
      scmd = ['sort','-k3,3','-k5,5n','-k6,6n','-k4,4']
    if tempdir: scmd += ['-T',tempdir.rstrip('/')+'/']
    if self._gz:
      cmd1 = "gzip"
      if os.name == 'nt':
         sys.stderr.write("WARNING: Windows OS Detected. close_fds not available.")
         p1 = Popen(cmd1,stdout=self._fh,stdin=PIPE,shell=True)
         p2 = Popen(" ".join(scmd),stdout=p1.stdin,stdin=PIPE,shell=True)
      else:
         p1 = Popen(cmd1.split(),stdout=self._fh,stdin=PIPE,close_fds=True)
         p2 = Popen(scmd,stdout=p1.stdin,stdin=PIPE,close_fds=True)
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
