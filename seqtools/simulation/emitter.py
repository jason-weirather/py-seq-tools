"""These classes should produce simulation products"""
from seqtools.simulation.randomsource import RandomSource

class TranscriptomeEmitter:
  """Give it a transcriptome definition and a reference genome for it
     initialy give it uniform probability

  :param transcriptome: A transcriptome from which to produce transcripts
  :param seed: Seeded random generation
  :param rand: A class that can generate randome numbers if you have one already seeded or want totally random
  :type transcriptome: Transcriptome
  :type seed: int
  :type rand: RandomSource
  """
  def __init__(self,transcriptome,seed=None,rand=None):
    if rand: self.random = rand
    elif seed: self.random = RandomSource(seed)
    else: self.random = RandomSource()

    self._transcriptome = transcriptome
    ######
    tcnt = len(self._transcriptome.get_transcripts())
    self._weights = [float(i+1)/float(tcnt) for i in range(0,tcnt)]
    ## _log stores what we are emitting ##
    self._log = []

  def emit_transcript(self):
    """Get a transcript according according to weight of transcript

    :return: One random Transcript
    :rtype: Transcript
    """
    i = self.random.get_weighted_random_index(self._weights)
    return self._transcriptome.get_transcripts()[i]

  def set_weights_by_dict(self,weights):
    """input: an array of weights <<txname1> <weight1>> <<txname2> <weight2>>...
       if this does not get set then even weighting will be used

    :param weights: [[tx1,wght1],[tx2,wght2],...[txN,wightN]]
    :type weights: list
    """
    self._weights = []
    txnames = [x.get_transcript_name() for x in self._transcriptome.get_transcripts()]
    for txname in txnames:
      if txname in weights:
        self._weights.append(float(weights[txname]))
      else:
        self._weights.append(float(0))
    return
