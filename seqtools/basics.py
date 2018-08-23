"""Put generally useful things here"""
import re, sys
from io import StringIO

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

class Capturing(list):
   """Capture stdout during part of code execution and store it like

   with Capturing() as output:
      do_stuff_that_goes_to_stdout()

   From   

   http://stackoverflow.com/questions/16571150/how-to-capture-stdout-output-from-a-python-function-call
   """
   def __enter__(self):
      self._cache = sys.stdout
      self._output = StringIO()
      sys.stdout = self._output
      return self
   def __exit__(self,*args):
      self.extend(self._output.getvalue().splitlines())
      del self._output
      sys.stdout = self._cache
