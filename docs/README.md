## AbiPy documentation

This is the top level build directory for the AbiPy documentation.  
All of the documentation is written using sphinx, a python documentation system built on top of ReST.  
This directory contains:

  * api - placeholders to automatically generate the api documentation

  * devel - documentation for AbiPy developers

  * sphinxext - Sphinx extensions for the AbiPy docs

  * users - the user documentation, e.g plotting tutorials, configuration tips, etc.

  * faq - frequently asked questions

  * index.rst - the top level include document for AbiPy docs

  * conf.py - the sphinx configuration

  * _static - used by the sphinx build system

  * _templates - used by the sphinx build system

To build the HTML documentation, install sphinx then type `make html` that will execute::

    sphinx-build -b html -d _build/doctrees . _build/html

Remeber to issue::

    export GENERATE_SPHINX_GALLERY=1

to activate the generation of the thumbnails in examples/flows.

The documentation is produced in `_build/html`.

You can run ``make help`` to see information on all possible make targets.

To install dependencies::

    pip install -R requirements.txt

To deploy to gh-pages::

   ./ghp_import.py _build/html/ -n -p
