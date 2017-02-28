"""These classes should produce simulation products"""
from seqtools.simulation.randomsource import RandomSource
from seqtools.simulation.permute import CutMaker
from seqtools.simulation.permute import ErrorMakerFlatRate
from collections import namedtuple
from seqtools.format.fastq import FASTQ

EmitterGenericOptions = namedtuple('ReadEmitterOptions',
   ['seed',
    'rand'])

class EmitterGeneric(object):
   """Emitters will vary based on the source they draw from"""
   def __init__(self,options=None):
      if not options: options = EmitterGeneric.Options()
      self._options = options
      #if self._options.rand: self._random = self._options.rand
      if not self._options.rand and self._options.seed: 
         self._options = self._options._replace(rand=RandomSource(self._options.seed))
      if not self._options.rand: 
         self._options = self._options._replace(rand=RandomSource())
      

   @property
   def options(self): return self._options

   def emit(self):
      raise ValueError('the class you are emitting from must override this emit')

   @staticmethod
   def Options(**kwargs):
      """ A method for declaring options for the class"""
      construct = EmitterGenericOptions #IMPORTANT!  Set this
      names = construct._fields
      d = {}
      for name in names: d[name] = None #default values
      for k,v in kwargs.iteritems():
         if k in names: d[k] = v
         else: raise ValueError('Error '+k+' is not a property of these options')
      """Create a set of options based on the inputs"""
      return construct(**d)
   

TranscriptomeEmitterOptions = namedtuple('TranscriptomeEmitterOptions',
   ['seed',
    'rand'])

class TranscriptomeEmitter(EmitterGeneric):
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
    super(TranscriptomeEmitter,self).__init__(options)
    #self._options = options

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

  #def emit_transcript(self):
  #  """Get a transcript according according to weight of transcript

  #  :return: One random Transcript
  #  :rtype: Transcript
  #  """
  #  i = self._random.get_weighted_random_index(self._weights)
  #  return self._transcriptome.transcripts[i]

  def emit(self):
    """Get a mapping from a transcript

    :return: One random Transcript sequence
    :rtype: sequence
    """
    i = self.options.rand.get_weighted_random_index(self._weights)
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
class ReadEmitter(EmitterGeneric):
   """Emit reads"""
   def __init__(self,emitter,options=None):
      if not options: options = ReadEmitter.Options()
      # get seed and rand from emitter if it wasnt specified
      if not options.seed and not options.rand:
         options = options._replace(seed=emitter.options.seed)
         options = options._replace(rand=emitter.options.rand)
      super(ReadEmitter,self).__init__(options)
      self._source = emitter
      self._cutter = CutMaker(rand=self.options.rand)
      self._errors = []
      #ErrorMakerFlatRate(rand=self.options.rand)
      self._flip = True

   @property
   def cutter(self): return self._cutter

   @property
   def errors(self): return self._errors

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

   def set_cut_maker(self,cutter):
      """this method sets a cutter to process the original 
         sequence into fragments"""
      self._cutter = cutter

   def add_error_maker(self,error_generator):
      """this method sets an object that will permute the sequence 
          with errors"""
      self._errors.append(error_generator)
   def reset_error_makers(self): self._errors = []

   def emit(self,rlen=150):
      """Emit a read based on a source sequence"""
      source_tx = self._source.emit()
      source_read = self._cutter.cut(source_tx)
      if self._flip and self.options.rand.random() < 0.5: source_read = source_read.rc()
      srname = self.options.rand.uuid4()
      seqfull = FASTQ('@'+self.options.rand.uuid4()+"\tlong\n"+str(source_read.sequence)+"\n+\n"+'I'*source_read.sequence.length+"\n")
      seqperm1 = seqfull.copy()
      seqperm2 = seqfull.copy()
      for e in self.errors:
        seqperm1 = e.permute(seqperm1)
        seqperm2 = e.permute(seqperm2)
      sleft = seqperm1[0:rlen]
      sleft = FASTQ('@'+sleft.name+"\tleft\n"+sleft.sequence+"\n+\n"+sleft.qual+"\n")
      sright = seqperm2.rc()[0:rlen]
      sright = FASTQ('@'+sright.name+"\tright\n"+sright.sequence+"\n+\n"+sright.qual+"\n")
      emission = TranscriptEmission(source_tx,
             Source(source_read,
                    source_read.slice_sequence(0,rlen),
                    source_read.rc().slice_sequence(0,rlen)),
             Read(seqperm1,
                  sleft,
                  sright
                 ))
      return emission

TranscriptEmission = namedtuple('TranscriptEmission',
   ['transcript',
    'source',
    'reads'])

Source = namedtuple('Source',
   ['fragment',
    'left',
    'right'])

Read = namedtuple('Read',
   ['fragment',
    'left',
    'right'])
