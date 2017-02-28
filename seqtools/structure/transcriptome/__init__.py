"""a module for working with families of transcripts

"""
import sys, random, string, uuid
from collections import namedtuple

class Transcriptome:
  """a class to store a transcriptome

  :param transcript_source: iterable transcript source
  :param ref_data:
  :type transcript_source: iterable 
  :type ref_data: dict()
  """
  def __init__(self,transcript_source=None,ref_fasta=None):
    self._transcripts = [x for x in transcript_source]
    self._ref = ref_fasta
    if ref_fasta:
      for i in range(0,len(self.transcripts)):
        self.transcripts[i].set_reference(self._ref)
        self.transcripts[i].sequence

  @property
  def transcripts(self):
    return self._transcripts
      
  def add_transcript(self,transcript):
    self.transcripts.append(transcript)

  def __str__(self):
    ostr = ''
    ostr += "Transcriptome containing "+str(len(self.transcripts))+" transcripts "
    ostr += "covering "+str(sum([x.length for x in self.transcripts]))+" bases"
    return ostr
