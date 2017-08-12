from setuptools import setup, find_packages
from codecs import open
from os import path
import sys

if sys.version_info < (2,7):
   sys.exit("Error: You are using Python "+str(sys.version_info)+"; Python 2.6 and below are not supported. Please use 2.7 or better\n")

this_folder = path.abspath(path.dirname(__file__))

setup(
  name='seq-tools',
  version='1.0.1',
  description='Python tools for working with biological sequence data',
  long_description='''
# py-seq-tools
Free python libraries for working with biological sequences, alignments, and formats. These libraries are part of Au-public scripts, I haven't documented them well yet.  The goal here is to put them into a package and get them better documented.

Longer term, these were written in python 2, but it would be nice to either support python3 full compatability or to completely migrate them over to python3.

The command line utilities of this package serve a dual purpose. 1) is to provide useful utilities you can run through the commandline, and 2) is to provide examples on how to make use the of the libraries included in this package.

Access to the command line utilities can be found with

`$ seq-tools <cli command>`

Where cli command is the name of a cli module.

To use the libraries in this distribution please consult the [API documentation](https://github.com/jason-weirather/py-seq-tools/blob/master/manual.pdf)

You can install seq-tools with pip

From within the source directory you can try

`$ pip install .`

This should make the `seq-tools` command available as well as the `seqtools` library for import.

This library includes the following modules:

* seqtools
  * align
  * errors
  * format
    * bamindex
    * bed
    * bgzf
    * fasta
    * fastq
    * pacbio
    * psl
    * sam
  * graph
  * range
  * sequence
  * simulation
    * emitter
    * permute
    * randomsource
  * statistics
  * stream
  * structure

Again, details are available in the [API documentation](https://github.com/jason-weirather/py-seq-tools/blob/master/manual.pdf).
  ''',
  url='https://github.com/jason-weirather/py-seq-tools',
  author='Jason L Weirather',
  author_email='jason.weirather@gmail.com',
  license='Apache License, Version 2.0',
  classifiers=[
    'Development Status :: 3 - Alpha',
    'Intended Audience :: Science/Research',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'License :: OSI Approved :: Apache Software License'
  ],
  keywords='bioinformatics, sequence, alignment',
  packages=['seqtools'],
  entry_points = {
    'console_scripts':['seq-tools=seqtools.cli.cli_front:main']
  }
)
