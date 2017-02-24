"""The testing suite for seqtools"""
import unittest, sys, gzip, hashlib, os
from tempfile import NamedTemporaryFile

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = THIS_DIR+'/../../data'

import seqtools.simulation
from seqtools.simulation.randomsource import RandomSource
class SimulationBasics(unittest.TestCase):
   """Test out various basic simulation functions"""
   def setUp(self):
      return
   def test_randomNT(self):
      """Test the seeded random"""
      r = RandomSource(1000)
      s = ''.join([r.random_nt() for x in range(1,10)])
      self.assertEqual('TGACCGTAG',s)
      s = ''.join([r.random_nt() for x in range(1,10)])
      self.assertEqual('CCAGACATC',s)
   def test_weightedRandomIndex(self):
      """Test the weighted random is working correctly"""
      v = [0.001,0.25,0.25,0.5,0.00001]
      r = RandomSource(200)
      cnt = 0
      tot = 20000
      for i in range(0,tot):
         if 3 == r.get_weighted_random_index(v): cnt+=1
      self.assertAlmostEqual(float(cnt)/tot,0.5,places=1)

from seqtools.structure.transcriptome import Transcriptome
from seqtools.format.gpd import GPDStream
from seqtools.format.fasta import FASTAData
from seqtools.simulation.emitter import TranscriptomeEmitter
from seqtools.simulation.emitter import ReadEmitter
class TryTranscriptomeEmitter(unittest.TestCase):
   def setUp(self):
      self.gpd_path = DATA_DIR+'/chr21chr22chrM.gencode25.gpd.gz'
      self.fasta_path = DATA_DIR+'/chr21chr22chrM.fa.gz'
      self.r = RandomSource(101)
      self.txome = Transcriptome(GPDStream(gzip.open(self.gpd_path)),
                                 FASTAData(gzip.open(self.fasta_path).read()))
      self.txe = TranscriptomeEmitter(self.txome,
                                 TranscriptomeEmitter.Options(rand=self.r))
   def test_txemitter(self):
      """Try to emit transcript sequences"""
      s = self.txe.emit_transcript().sequence
      
      h = hashlib.md5(str(s)).hexdigest()

      self.assertEqual(h,'5dfe012a446262d7a7703937cd340ebf')
      s = self.txe.emit_transcript().sequence
      h = hashlib.md5(str(s)).hexdigest()
      self.assertEqual(h,'c902ff6bf269afecf5411801f99dba36')
   def test_reademitter(self):
      """Test emitting reads"""
      re = ReadEmitter(self.txome)

if __name__ == '__main__':
   unittest.main(verbosity=2)
