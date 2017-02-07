"""Classes to work with sam headers"""
import sys

class SAMHeader:
  """class to retain information about a SAMheader and access data from it"""
  def __init__(self,header_text):
    self._text = header_text
    self.tags = []
    self._sequence_lengths = {}
    for line in self._text.split("\n"):
      if len(line) == 0: continue
      tag = line[1:3]
      rem = line[4:]
      self.tags.append({'tag':tag,'info':{}})
      for c in [{'field':x[0:2],'value':x[3:]} for x in rem.split("\t")]:
        self.tags[-1]['info'][c['field']] = c['value']
    for v in [x['info'] for x in self.tags if x['tag'] == 'SQ']:
      self._sequence_lengths[v['SN']] = int(v['LN'])
    return

  def __str__(self): return self._text.rstrip()

  def text(self): return self._text.rstrip()+"\n"

  @property
  def sequence_names(self):
    return self._sequence_lengths.keys()

  @property
  def sequence_lengths(self):
    """return a dictionary of sequence lengths"""
    return self._sequence_lengths

  def get_sequence_length(self,sname):
    """get a specific sequence length"""
    return self._sequence_lengths[sname]
