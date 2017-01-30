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
