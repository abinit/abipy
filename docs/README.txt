abipy documentation
===================

This is the top level build directory for the abipy documentation.  
All of the documentation is written using sphinx, a
python documentation system built on top of ReST.  
This directory contains:

  * api - placeholders to automatically generate the api documentation

  * devel - documentation for abipy developers

  * sphinxext - Sphinx extensions for the abipy docs

  * users - the user documentation, e.g plotting tutorials, configuration tips, etc.

  * faq - frequently asked questions

  * index.rst - the top level include document for abipy docs

  * conf.py - the sphinx configuration

  * _static - used by the sphinx build system

  * _templates - used by the sphinx build system

To build the HTML documentation, install sphinx (v1.1.3 or greater required), 
then type "make html" that will execute: 

  sphinx-build -b html -d _build/doctrees . _build/html

The documentation is produced in _build/html.
