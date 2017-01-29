from setuptools import setup, find_packages
from codecs import open
from os import path

this_folder = path.abspath(path.dirname(__file__))
with open(path.join(this_folder,'README.md'),encoding='utf-8') as inf:
  long_description = inf.read()

setup(
  name='seq-tools',
  version='0.2.0',
  description='Python tools for working with biological sequence data',
  long_description=long_description,
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
  #packages=find_packages('seqtools'),
  packages=['seqtools'],
  #package_directories={'',''},
  entry_points = {
    'console_scripts':['seq-tools=seqtools.cli.cli_front:main']
  }
)
  
