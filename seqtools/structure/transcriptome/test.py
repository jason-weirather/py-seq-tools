"""The testing suite for seqtools"""
import unittest, sys, gzip, hashlib, os
from tempfile import NamedTemporaryFile

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = THIS_DIR+'/../../../data'
from seqtools.format.gpd import GPDStream
from seqtools.format.fasta import FASTAData
import seqtools.structure.transcriptome
class Transcriptome(unittest.TestCase):
   """Test out various simulation functions"""
   def setUp(self):
      self.gpd_path = DATA_DIR+'/chr21chr22chrM.gencode25.gpd.gz'
      self.fasta_path = DATA_DIR+'/chr21chr22chrM.fa.gz'
   def test_makeTranscriptome(self):
      """Try loading a transcriptome with a fasta and genepred"""
      gpds = GPDStream(gzip.open(self.gpd_path))
      fasta = FASTAData(gzip.open(self.fasta_path).read())
      txome = seqtools.structure.transcriptome.Transcriptome(gpds,fasta)
      neg_tx = [x for x in txome.transcripts if x.strand=='-']
      pos_tx = [x for x in txome.transcripts if x.strand=='+']
      last_tx = txome.transcripts[-1]
      self.assertEqual(len(neg_tx),3360)
      self.assertEqual(len(pos_tx),3562)
      hash =  hashlib.md5(str(last_tx.sequence).encode('utf-8')).hexdigest()
      val2 = '9b7cf1954606443259d8ae3fc4063e4d'
      self.assertEqual(val2,hash)
if __name__ == '__main__':
   unittest.main(verbosity=2)
