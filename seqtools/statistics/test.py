"""The testing suite for seqtools"""
import unittest, sys, gzip, hashlib, cStringIO
from seqtools.range import GenomicRange
from tempfile import NamedTemporaryFile

import seqtools.statistics as statistics
class Statistics(unittest.TestCase):
   """ Start with basic tests of statistics"""
   def setUp(self):
      """Testing the basic statistics functions"""
      self.dat = [1,2,3,4,5,6,9,10,11,12,13,14]

   def test_median(self):
      """Simple median"""
      self.assertAlmostEqual(
         statistics.median(self.dat),
         7.5)

   def test_average(self):
      """Simple mean (average)"""
      self.assertAlmostEqual(
         statistics.average(self.dat),
         7.5)

   def test_stddev(self):
      """Simple standard deviation"""
      v = statistics.standard_deviation(self.dat)
      self.assertAlmostEqual(v,4.54272645405)

   def test_N50(self):
      """N50 is used in genome assembly assessment"""
      v = statistics.N50(self.dat)
      self.assertEqual(v,11)

if __name__ == '__main__':
   unittest.main(verbosity=2)
