import re

class PacBioReadName:
  """A class to get information from PacBio read names
     Only supports subreads and ccs reads for now, but may be
     nice to add isoseq support in the future

  :param name: Name is a pacbio name
  :type name: string
  """
  def __init__(self,name):
    self._name = name.rstrip()
    m1 = re.match('([^\/]+)\/(\d+)\/ccs$',self._name)
    m2 = re.match('([^\/]+)\/(\d+)\/(\d+)_(\d+)$',self._name)
    self._base = None
    self._molecule_number = None
    if m1:
      self._is_ccs = True
      self._base = m1.group(1)
      self._molecule_number = int(m1.group(2))
    elif m2:
      self._is_ccs = False
      self._base = m2.group(1)
      self._molecule_number = int(m2.group(2))
      self._sub_start = int(m2.group(3))
      self._sub_end = int(m2.group(4))
    if not self._base: return None
  def __str__(self):
    """ The string representation of this is just the name """
    return self._name
  def name(self):
    """Just return the full name of the read"""
    return self._name;

  def is_ccs(self):
    """ Return true if this sequence is a ccs read"""
    if self._is_ccs: return True
    return False

  def is_sub(self):
    """ Return true if this sequence is a subread"""
    if self._is_sub: return True
    return False;

  def get_base(self):
    """get the first part of the name specific to the SMRT cell"""
    return self._base;

  def get_molecule(self):
    return self._base+'/'+str(self._molecule_number)

  def get_molecule_number(self):
    """get the second part of the name specific to the molecule in the flow cell."""
    return self._molecule_number
