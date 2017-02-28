"""A class to aid in generating random numbers and sequences

It doesn't seem necessary to create an options class since this class will probably not be extended

"""
import random, sys

nts = ['A','C','G','T']

hexchars = '0123456789abcdef'
uuid4special = '89ab'

class RandomSource:
  """You can asign it a seed if you want

  :param seed: seed the pseduorandom number generator
  :type seed: int
  """
  def __init__(self,seed=None):
    self._random = random.Random()
    if seed: self._random.seed(seed)

  def choice(self,arr):
    """Uniform random selection of a member of an list

    :param arr: list you want to select an element from
    :type arr: list
    :return: one element from the list
    """
    ind = self.randint(0,len(arr)-1)
    return arr[ind]
  
  def random(self):
    """generate a random number

    :return: uniform random float between 0 and 1
    :rtype: float
    """
    return self._random.random()

  def gauss(self,mu,sigma):
    """Generate a random number based on a gaussian distribution

    :param mu: mean of distribution
    :param sigma: standard deveiation of distribution (i think)
    :type mu: float
    :type sigma: float
    """
    return self._random.gauss(mu,sigma)

  def randint(self,a,b):
    """Generate a random integer uniform distribution between a and b like randint of the usual random class

    :return: random int between a and b
    :rtype: int
    """
    return self._random.randint(a,b)

  def different_random_nt(self,nt):
    global nts
    """generate a random nucleotide change. uniform random.  will never return itself

    :param nt: current nucleotide
    :type nt: char
    :return: new nucleotide
    :rtype: char
    """
    return self._random.choice([x for x in nts if x != nt.upper()])

  def random_nt(self):
    """Produce a random nucleotide (uniform random)

    :return: nucleotide
    :rtype: char
    """
    global nts
    return self._random.choice(nts)
  
  def get_weighted_random_index(self,weights):
    """Return an index of an array based on the weights
       if a random number between 0 and 1 is less than an index return the lowest index

    :param weights: a list of floats for how to weight each index [w1, w2, ... wN]
    :type weights: list
    :return: index
    :rtype: int
    """
    tot = float(sum([float(x) for x in weights]))
    fracarray = [weights[0]]
    for w in weights[1:]:
      prev = fracarray[-1]
      fracarray.append(w+prev)
    #print fracarray
    rnum = self._random.random()*tot
    #print rnum
    #sys.exit()
    for i in range(len(weights)):
      if rnum < fracarray[i]: return i
    sys.stderr.write("Warning unexpected no random\n")

  def uuid4(self):
    """Make an id in the format of UUID4, but keep in mind this could very well be pseudorandom, and if it is you'll not be truely random, and can regenerate same id if same seed"""
    return ''.join([hexchars[self.randint(0,15)] for x in range(0,8)]) + '-' +\
           ''.join([hexchars[self.randint(0,15)] for x in range(0,4)]) + '-' +\
           '4'+''.join([hexchars[self.randint(0,15)] for x in range(0,3)]) + '-' +\
           uuid4special[self.randint(0,3)]+''.join([hexchars[self.randint(0,15)] for x in range(0,3)]) + '-' +\
           ''.join([hexchars[self.randint(0,15)] for x in range(0,12)])
    
