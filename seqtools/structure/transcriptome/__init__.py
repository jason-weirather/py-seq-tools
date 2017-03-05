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
    self._transcripts = []
    if transcript_source:
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

  def sort_transcripts(self):
     """Sort the transcripts stored here"""
     txs = sorted(self.transcripts,key=lambda x: (x.range.chr, x.range.start, x.range.end))
     self._transcripts = txs

  def __str__(self):
    ostr = ''
    ostr += "Transcriptome containing "+str(len(self.transcripts))+" transcripts "
    ostr += "covering "+str(sum([x.length for x in self.transcripts]))+" bases"
    return ostr

  def transcript_stream(self):
     for tx in self._transcripts:
        yield tx
     
