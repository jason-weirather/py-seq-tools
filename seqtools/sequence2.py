import re, sys, base64, zlib
from string import maketrans

class Seq:
  """Basic Sequence structure

  :param seq: nucleotide sequence
  :param name: name is optional
  :type seq: string
  :type name:

  """
  def __init__(self,seq=None,name=None):
    self.name = name
    self.seq = seq
  def __getitem__(self,key):
    #if its a slice deal with the sequence only
    if isinstance(key,slice):
      newseq = self.seq[key.start:min(key.stop,len(self.seq))]
      return Seq(newseq,self.name)
    return {'name':self.name,'seq':self.seq}[key]
  def __str__(self):
    return self.seq
  def __len__(self):
    return len(self.seq)
  def rc(self):
    """reverse complement

    :return: reverse complemented sequence of same name
    :rtype: Seq
    """
    return Seq(rc(self.seq),self.name)
  def copy(self):
    """a new sequence object

    :return: copy of the sequence
    :rtype: Seq

    """
    return Seq(self.seq,self.name)


  def fasta(self):
    """Get the fasta formated string
    Pre: seq and name are set
    Post: string representation of the fasta entry

    :return: fasta string
    :rtype: string

    """
    return '>'+self.name+"\n"+self.seq+"\n"

  def gc_content(self):
    """Calculate the GC content of the sequence

    :return: GC content
    :rtype: float

    """
    if len(self.seq) == 0: return None
    n_count = self.n_count()
    if len(self.seq) - n_count == 0: return None
    return float(self.seq.translate(maketrans('GCgc','GGGG')).count('G')-n_count)/float(len(self.seq)-n_count)

  def n_count(self):
    """ Count the numbers of 'N's in a sequence.  case insensitive.

    :return: N count
    :rtype: int
    """
    return self.seq.translate(maketrans('Nn','NN')).count('N')

def rc(seq):
  """Fast reverse complement function using a translation table and slice

  :param seq: string to reverse complement
  :type seq: string
  :return: reverse complemented sequence
  :rtype: string
  """
  complement = maketrans('ACTGUNXactgunx','TGACANXtgacanx')
  return seq.translate(complement)[::-1]


def encode_name(conversion_string):
  """Make a name into an encoding that can store any character

  :param conversion_string: thing to be encoded
  :type conversion_string: string
  :return: encoded_name
  :rtype: string
  """
  compressed_string = zlib.compress(conversion_string,9)
  enc_string = base64.b32encode(compressed_string)
  return 'SZ_'+enc_string.rstrip('=')

def decode_name(safename):
  """Make an encoded name into its decoded.

  :param safename: thing to be decoded
  :type safenameg: string
  :return: decoded name
  :rtype: string
  """
  frag = safename.lstrip('SZ_')
  padding = (8-(len(frag) % 8)) % 8
  c = base64.b32decode(frag+'='*padding)
  return  zlib.decompress(c)

