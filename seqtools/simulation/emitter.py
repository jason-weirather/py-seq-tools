"""These classes should produce simulation products"""
from seqtools.simulation.randomsource import RandomSource
from collections import namedtuple

TranscriptomeEmitterOptions = namedtuple('TranscriptomeEmitterOptions',
   ['seed',
    'rand'])

class TranscriptomeEmitter(object):
  """Give it a transcriptome definition and a reference genome for it
     initialy give it uniform probability

  :param transcriptome: A transcriptome from which to produce transcripts
  :param seed: Seeded random generation
  :param rand: A class that can generate randome numbers if you have one already seeded or want totally random
  :type transcriptome: Transcriptome
  :type seed: int
  :type rand: RandomSource
  """
  def __init__(self,transcriptome,options=None):
    if not options: options = TranscriptomeEmitter.Options()
    self._options = options
    if self._options.rand: self._random = self._options.rand
    elif self._options.seed: self._random = RandomSource(seed)
    else: self._random = RandomSource()

    self._transcriptome = transcriptome
    ######
    tcnt = len(self._transcriptome.transcripts)
    self._weights = [float(i+1)/float(tcnt) for i in range(0,tcnt)]
    ## _log stores what we are emitting ##
    self._log = []

  @staticmethod
  def Options(**kwargs):
     """ A method for declaring options for the class"""
     construct = TranscriptomeEmitterOptions #IMPORTANT!  Set this
     names = construct._fields
     d = {}
     for name in names: d[name] = None #default values
     for k,v in kwargs.iteritems():
       if k in names: d[k] = v
       else: raise ValueError('Error '+k+' is not a property of these options')
     """Create a set of options based on the inputs"""
     return construct(**d)

  def emit_transcript(self):
    """Get a transcript according according to weight of transcript

    :return: One random Transcript
    :rtype: Transcript
    """
    i = self._random.get_weighted_random_index(self._weights)
    return self._transcriptome.transcripts[i]

  def set_weights_by_dict(self,weights):
    """input: an array of weights <<txname1> <weight1>> <<txname2> <weight2>>...
       if this does not get set then even weighting will be used

    :param weights: [[tx1,wght1],[tx2,wght2],...[txN,wightN]]
    :type weights: list
    """
    self._weights = []
    txnames = [x.get_transcript_name() for x in self._transcriptome.transcripts]
    for txname in txnames:
      if txname in weights:
        self._weights.append(float(weights[txname]))
      else:
        self._weights.append(float(0))
    return

ReadEmitterOptions = namedtuple('ReadEmitterOptions',
   ['seed',
    'rand'])

class ReadEmitter(TranscriptomeEmitter):
   """Emit reads"""
   def __init__(self,tx_emitter,options=None):
      if not options: options = ReadEmitter.Options()
      super(ReadEmitter,self).__init__(tx_emitter,options)

   @staticmethod
   def Options(**kwargs):
      """ A method for declaring options for the class"""
      construct = ReadEmitterOptions #IMPORTANT!  Set this
      names = construct._fields
      d = {}
      for name in names: d[name] = None #default values
      for k,v in kwargs.iteritems():
         if k in names: d[k] = v
         else: raise ValueError('Error '+k+' is not a property of these options')
      """Create a set of options based on the inputs"""
      return construct(**d)

   
