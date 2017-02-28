"""The testing suite for seqtools"""
import unittest, sys, gzip, hashlib, os
from tempfile import NamedTemporaryFile

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = THIS_DIR+'/../../data'
from seqtools.format.gpd import GPDStream
from seqtools.format.fasta import FASTAData
import seqtools.structure.transcriptome
class Structure(unittest.TestCase):
   """Test out various simulation functions"""
   def setUp(self):
      self.gpd_path = DATA_DIR+'/chr21chr22chrM.gencode25.gpd.gz'
      self.fasta_path = DATA_DIR+'/chr21chr22chrM.fa.gz'
      """Try loading a transcriptome with a fasta and genepred"""
      self.gpds = list(GPDStream(gzip.open(self.gpd_path)))
      self.fasta = FASTAData(gzip.open(self.fasta_path).read())
      #txome = seqtools.structure.transcriptome.Transcriptome(gpds,fasta)
   def test_sequencelength(self):
      l =  self.gpds[2].length
      self.assertEqual(l,3412)
   def test_sequence(self):
      self.gpds[2].set_reference(self.fasta)
      s = self.gpds[2].sequence
      self.assertEqual(hashlib.md5(str(s)).hexdigest(),'7d90fcf2afa796e5c6eebd10f1cf6433')
   def test_targetslice(self):
      gpd = self.gpds[2]
      g2 = gpd.slice_target('chr21',26503701,26567958) # save the last base of the second exon
      self.assertEqual(g2.exons[0].length,1) # one base left in first exon
      self.assertEqual(g2.exons[-1].length,2) #two bases left in last exon
      self.assertEqual(len(g2.exons),3)

   def test_sequenceslice(self):
      gpd = self.gpds[2]
      gpd.set_reference(self.fasta)
      g2 = gpd.slice_sequence(182,3000) # save the last base of the second exon
      self.assertEqual(str(g2.sequence),str(gpd.sequence[182:3000]))

if __name__ == '__main__':
   unittest.main(verbosity=2)
