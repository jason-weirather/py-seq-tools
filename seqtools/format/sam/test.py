"""The testing suite for seqtools"""
import unittest, sys, gzip, hashlib, io, os
from seqtools.range import GenomicRange
from tempfile import NamedTemporaryFile
from seqtools.format.fasta import FASTAData

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = THIS_DIR+'/../../../data'

import seqtools.cli.utilities.sort as cmd_sort
class SortCmd(unittest.TestCase):
   def test_samsortcmd(self):
      """Try and sort a sam file"""
      rawpath = DATA_DIR+'/chr21chr22chrM.bam'
      f = NamedTemporaryFile(delete=True,suffix='.sam')
      cmd_sort.external_cmd('sort '+rawpath+' --bam --threads 2 -o '+f.name)
      hash = hashlib.md5(open(f.name).read()).hexdigest()
      self.assertEqual('e9740957cea11c6e338a1edfdbe51d68',hash)

import seqtools.format.sam as sam
import seqtools.format.sam.bam.files as bamfiles
class ALIGN(unittest.TestCase):
   """Test out various sam functions"""
   def setUp(self):
      path_ref = DATA_DIR+'/chr21chr22chrM.fa.gz'
      data_ref = gzip.open(path_ref).read()
      self.ref = FASTAData(data_ref)
      self.rawpath = DATA_DIR+'/chr21chr22chrM.bam'
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
      stream = io.StringIO(self.sam)
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

if __name__ == '__main__':
   unittest.main(verbosity=2)
