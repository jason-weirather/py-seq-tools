"""Put generally useful things here"""
import re

def is_uuid4(instr):
  """A validator to confirm a string is indeed a UUID4

  :param instr: input string that may be uuid4
  :type instr: String
  :return: true if it is in UUID4 format
  :rtype: bool
  """
  v = instr.strip().replace('-','').lower()
  if len(v) != 32: return False
  if not re.match('^[0-9a-f]+$',v): return False
  if v[12] != '4': return False
  if not re.match('[89ab]',v[16]): return False
  return True
