"""The testing suite for seqtools"""
import unittest, sys, gzip, hashlib, StringIO
from seqtools.range import GenomicRange

import seqtools.format.sam as sam
import seqtools.format.sam.bam.files as bamfiles
class ALIGN(unittest.TestCase):
   """Test out various sam functions"""
   def setUp(self):
      path_ref = '../data/chr21chr22chrM.fa.gz'
      data_ref = gzip.open(path_ref).read()
      self.ref = fasta.FASTAData(data_ref)
      self.rawpath = '../data/chr21chr22chrM.bam'
      bf = bamfiles.BAMFile(self.rawpath)
      self.sam = bf.header.text.rstrip()+"\n"
      for b in bf:  
         self.sam += b.sam_line+"\n"
      self.assertEqual(16569,bf.header.sequence_lengths['chrM'])
   def test_bamread(self):
      """Test the bam was read in properly and converted to sam properly"""
      samhash = hashlib.md5(self.sam).hexdigest()
      self.assertEqual('2e92a4793a4a2140103dcb26d2523d68',samhash)
   def test_samread(self):
      """Test that the sam can be streamed properly"""
      stream = StringIO.StringIO(self.sam)
      ss = sam.SAMStream(stream)
      buffer = ss.header.text.rstrip()+"\n"
      for s in ss:
         buffer += s.sam_line+"\n"
      samhash = hashlib.md5(buffer).hexdigest()
      self.assertEqual('2e92a4793a4a2140103dcb26d2523d68',samhash)
      self.assertEqual(16569,ss.header.sequence_lengths['chrM'])
   def test_pslconversion(self):
      """Check conversion to PSL"""
      bf = bamfiles.BAMFile(self.rawpath,bamfiles.BAMFile.Options(reference=self.ref))
      buffer = ''
      of = open('check','w') 
      for b in bf:
         if not b.is_aligned(): continue  
         p =  b.get_PSL()
         buffer += str(p)+"\n"
      self.assertEqual('56e954a5efd405ab1449c0af746846d6',hashlib.md5(buffer).hexdigest())

import seqtools.format.fasta as fasta
class FASTA(unittest.TestCase):
   def setUp(self):
      """Loading in the test data."""
      self.rawpath = '../data/chr21chr22chrM.fa.gz'
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
