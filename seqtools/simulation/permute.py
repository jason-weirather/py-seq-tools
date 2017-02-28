"""These classes are here to alter simulation products or other sequences"""
import math
from seqtools.sequence import rc, Sequence
from seqtools.simulation.randomsource import RandomSource
from seqtools.format.fastq import FASTQ

class ErrorMakerGeneric(object):
   """Add errors to sequences"""
   def __init__(self,rand=None,seed=None):
      if rand:
         self.random = rand
      else:
         if seed: self.random = RandomSource(seed)
         else: self.random = RandomSource()
   def permute(self,seq):
      return seq

class ErrorMakerFlatRate(ErrorMakerGeneric):
  """Class to define how to make errors, and to introduce those errors"""
  def __init__(self,rate=0,rand=None,seed=None):
     super(ErrorMakerFlatRate,self).__init__(rand=rand,seed=seed)
     self._rate = rate
  def set_rate(self,rate): self._rate = rate
  def permute(self,fastq):
     sequence = fastq.sequence
     seq = ''
     qual = ''
     qbase = rate_to_phred33(self._rate)
     for i in range(len(sequence)):
        # check context
        rnum = self.random.random()
        if rnum < self._rate:
           step1 = self._rate/float(5)
           if rnum < 3*step1:
              seq += self.random.different_random_nt(sequence[i])
              qual += qbase
           elif rnum < 4*step1:
              if (rnum-3*step1)/step1 < 0.5: #either before or after
                 seq = seq+sequence[i]+self.random.random_nt()
                 qual += qbase+qbase
              else:
                 seq = seq+self.random.random_nt()+sequence[i]
                 qual += qbase+qbase
           #we don't need to explicity do the deletion.
        else:
           seq += sequence[i]
           qual += qbase
     return FASTQ('@'+fastq.header+"\n"+seq+"\n+\n"+qual+"\n")
     
class ErrorMaker(ErrorMakerGeneric):
  """Class to define how to make errors, and to introduce those errors"""
  def __init__(self,rand=None,seed=None):
    super(ErrorMaker,self).__init__(rand=rand,seed=seed)

    #### context information ####
    self._before_base = None
    self._after_base = None
    #### set the reference base to change for del,mismatch ###
    self._observed_base = None
    #### set waht to change base to for ins or mismatch
    self._modified_base = None
    self._substitution_rate = 0
    self._deletion_rate = 0
    self._insertion_rate = 0

  def set_substitution_rate(self,rate):
    self._substitution_rate = rate
  def set_deletion_rate(self,rate):
    self._deletion_rate = rate
  def set_insertion_rate(self,rate):
    self._insertion_rate = rate

  def permute(self,fastq):
    s = fastq
    if self._substitution_rate > 0:
      s = self.random_substitution(s, self._substitution_rate)
    if self._insertion_rate > 0:
      s = self.random_insertion(s, self._insertion_rate)
    if self._deletion_rate > 0:
      s = self.random_deletion(s, self._deletion_rate)
    return s

  def set_before_context(self,base):
    """Limit errors to a specific preceeding base context"""
    self._before_base = base
  def set_after_context(self,base):
    """Limit errors to a specific following base context"""
    self._after_base = base
  def set_observed_base(self,base):
    """Limit errors to a specific reference base"""
    self._observed_base = base
  def set_modified_base(self,base):
    """Limit errors to a specific type of sequenced base"""
    self._modified_base = base

  def random_substitution(self,fastq,rate):
    """Perform the permutation on the sequence

    :param fastq: FASTQ sequence to permute
    :type fastq: format.fastq.FASTQ
    :param rate: how frequently to permute
    :type rate: float
    :return: Permutted FASTQ
    :rtype: format.fastq.FASTQ
    """
    sequence = fastq.sequence
    seq = ''
    for i in range(len(sequence)):
      # check context
      prev = None
      if i >= 1: prev = sequence[i-1]
      next = None
      if i < len(sequence)-1: next = sequence[i+1]
      if self._before_base and (not prev or prev != self._before_base): 
        seq+=sequence[i]
        continue
      if self._after_base and (not next or next != self._after_base): 
        seq+=sequence[i]
        continue
      if self._observed_base and (sequence[i] != self._observed_base):
        seq+=sequence[i]
        continue

      rnum = self.random.random()
      if rnum < rate:
        if not self._modified_base:
          seq += self.random.different_random_nt(sequence[i])
        else:
          seq += self._modified_base
      else:
        seq += sequence[i]
    return FASTQ('@'+fastq.header+"\n"+seq+"\n+\n"+fastq.qual+"\n")

  def random_deletion(self,fastq,rate):
    """Perform the permutation on the sequence

    :param fastq: FASTQ sequence to permute
    :type fastq: format.fastq.FASTQ
    :param rate: how frequently to permute
    :type rate: float
    :return: Permutted FASTQ
    :rtype: format.fastq.FASTQ
    """
    sequence = fastq.sequence
    quality = fastq.qual
    seq = ''
    qual = None
    if quality: qual = ''
    for i in range(len(sequence)):
      # check context
      prev = None
      if i >= 1: prev = sequence[i-1]
      next = None
      if i < len(sequence)-1: next = sequence[i+1]
      if self._before_base and (not prev or prev != self._before_base): 
        seq+=sequence[i]
        if quality: qual+=quality[i]
        continue
      if self._after_base and (not next or next != self._after_base): 
        seq+=sequence[i]
        if quality: qual+=quality[i]
        continue
      if self._observed_base and (sequence[i] != self._observed_base):
        seq+=sequence[i]
        if quality: qual+=quality[i]
        continue

      rnum = self.random.random()
      if rnum >= rate:
        seq += sequence[i]
        if quality: qual+=quality[i]
    return FASTQ('@'+fastq.header+"\n"+seq+"\n+\n"+qual+"\n")

  def random_insertion(self,fastq,rate,max_inserts=1):
    """Perform the permutation on the sequence. If authorized to do multiple bases they are done at hte rate defined here.

    :param fastq: FASTQ sequence to permute
    :type fastq: format.fastq.FASTQ
    :param rate: how frequently to permute
    :type rate: float
    :param max_inserts: the maximum number of bases to insert (default 1)
    :type rate: int
    :return: Permutted FASTQ
    :rtype: format.fastq.FASTQ
    """
    sequence = fastq.sequence
    quality = fastq.qual
    seq = ''
    qual = None
    ibase = rate_to_phred33(rate)
    if quality: qual = ''
    z = 0
    while self.random.random() < rate and z < max_inserts:
      if self._before_base: break # can't do this one
      if self._after_base:
        if self._after_base != sequence[1]: break
      z += 1
      if self._modified_base:
        seq += self._modified_base
        if quality: qual += ibase
      else:
        seq += self.random.random_nt()
        if quality: qual += ibase
    z = 0
    for i in range(len(sequence)):
      # check context
      prev = sequence[i]
      next = None
      if i < len(sequence)-1: next = sequence[i+1]
      if self._before_base and (not prev or prev != self._before_base): 
        seq+=sequence[i]
        if quality: qual+=quality[i]
        continue
      if self._after_base and (not next or next != self._after_base): 
        seq+=sequence[i]
        if quality: qual+= quality[i]
        continue

      seq += sequence[i]
      if quality: qual += quality[i]
      while self.random.random() < rate and z < max_inserts:
        z+=1
        if self._modified_base:
          seq += self._modified_base
          if quality: qual += ibase
        else:
          seq += self.random.random_nt()
          if quality: qual += ibase
      z = 0
    return FASTQ('@'+fastq.name+"\n"+seq+"\n+\n"+qual+"\n")

  def random_flip(self,sequence):
    """Change the direction of the sequence with 0.5 probability"""
    if self.random.random() < 0.5:
      return rc(sequence)
    return sequence

class CutMaker:
  """Class to cut the sequence to different sizes

  :param rand: pass a random source, otherwise it gets a new RandomSource
  :param seed: if you want to set a seed here
  :type rand: RandomSource
  :type seed: int
  """
  def __init__(self,rand=None,seed=None):
    #print rand
    if rand:
      self.random = rand
    else:
      self.random = RandomSource()
      if seed: self.random = RandomSource(seed)
    self._gauss_min = None
    self._gauss_mu = None
    self._gauss_sigma = None
    #self.set_lr_cuts()

  def cut(self,seq):
    #cut a sequence or a mapping
    if not self._gauss_min: return seq # if not set up return the sequence
    rgauss = self.random.gauss(self._gauss_mu,self._gauss_sigma)
    l = min(seq.length,max(self._gauss_min,int(rgauss)))
    #print self._gauss_min
    #print self._gauss_mu
    #print rgauss
    leeway = seq.length-l
    start = self.random.randint(0,leeway)
    return seq.slice_sequence(start,start+l)

  def set_custom(self,gmin,gmu,gsigma):
    """Set a minimum lengtha, and then the gaussian distribution parameters for cutting
       For any sequence longer than the minimum the guassian parameters will be used"""
    self._gauss_min = gmin
    self._gauss_mu = gmu
    self._gauss_sigma = gsigma
  def set_lr_cuts(self):
    self._gauss_min = 1000
    self._gauss_mu = 4000
    self._gauss_sigma = 500
  def set_sr_cuts(self):
    self._gauss_min = 150
    self._gauss_mu = 290
    self._gauss_sigma = 290

def random_flip(sequence,rnum=None):
  """Flip a sequence direction with 0.5 probability"""
  randin = rnum
  if not randin: randin = RandomSource()
  if randin.random() < 0.5:
    return rc(sequence)
  return sequence

def rate_to_phred33(rate):
  """Convert an error rate to a phred 33 character"""
  if rate < 0.0001: return 'I'
  return chr(int(-10*math.log10(rate))+33)

def phred33_to_rate(q):
  """Convert a phred33 character to an error rate"""
  return math.pow(10,float(ord(q)-33)/-10)
