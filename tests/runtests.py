"""The testing suite for seqtools

Testing modules are amongst the module in test.py files.
Each of these modules points to the base data directory as needed

This script is just for executing the test run

"""
import unittest, sys
import seqtools.format.sam.test
import seqtools.format.fasta.test
import seqtools.statistics.test
import seqtools.simulation.test
import seqtools.structure.transcriptome.test
import seqtools.structure.test

if __name__ == '__main__':
   loader = unittest.TestLoader()
   s = []
   #s.append(loader.loadTestsFromModule(seqtools.format.fasta.test)) 
   #s.append(loader.loadTestsFromModule(seqtools.statistics.test)) 
   #s.append(loader.loadTestsFromModule(seqtools.format.sam.test)) 
   #s.append(loader.loadTestsFromModule(seqtools.structure.transcriptome.test)) 
   #s.append(loader.loadTestsFromModule(seqtools.structure.test)) 
   s.append(loader.loadTestsFromModule(seqtools.simulation.test)) 
   unittest.TextTestRunner(verbosity=2).run(unittest.TestSuite(s))
