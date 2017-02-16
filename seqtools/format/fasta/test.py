"""The testing suite for seqtools"""
import unittest, sys, gzip, hashlib, cStringIO, os
from seqtools.range import GenomicRange
from tempfile import NamedTemporaryFile

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = THIS_DIR+'/../../../data'

import seqtools.format.fasta as fasta
class FASTA(unittest.TestCase):
   def setUp(self):
      """Loading in the test data."""
      self.rawpath = DATA_DIR+'/chr21chr22chrM.fa.gz'
      self.rawdata = gzip.open(self.rawpath).read()
      self.ref = fasta.FASTAData(self.rawdata)
   def test_rawdata(self):
      """Make sure the testing data is what we expect"""
      v = 'cef4059e8384565ad437a4fb6a883045'
      self.assertEqual(hashlib.md5(self.rawdata).hexdigest(),v)
   def test_keys(self):
      """Test access by keys FASTAData"""
      out ='chr21,chr22,chrM'
      chrs = ','.join([x for x in self.ref])
      self.assertEqual(chrs,out) # access by iterator
      self.assertEqual(','.join(self.ref.keys()),out) #access by keys
   def test_length(self):
      """Test length of sequence FASTAData"""
      self.assertEqual(self.ref['chrM'].length,16569)
      self.assertEqual(len(self.ref['chrM']),16569)
      self.assertEqual(len(str(self.ref['chrM'])),16569)
   def test_preservation(self):
      """Make sure format is preserved"""
      v = ">myval\tcheck 1\tcheck2\nAAttCCgg\nATCG\nAAATTTGNnT\n"
      vstring = 'AAttCCggATCGAAATTTGNnT'
      fv = fasta.FASTA(v)
      self.assertEqual(fv.FASTA(),v)
      self.assertEqual(str(fv),vstring)
      self.assertEqual(fv.header,"myval\tcheck 1\tcheck2")
      self.assertEqual(fv.name,"myval")
   def test_streaming(self):
      """Make sure our streaming iterator is working"""
      inf = gzip.open(self.rawpath)
      stream = fasta.FASTAStream(inf)
      chrs = []
      for val in stream:
         chrs.append(val.name)
      self.assertEqual('chrM',chrs[-1])

if __name__ == '__main__':
   unittest.main(verbosity=2)
