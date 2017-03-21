#!/bin/bash
set -ev  # exit on first error, print each command
conda install sphinx
#conda install sphinx_rtd_theme
#conda install ipython
pip install sphinxcontrib.napoleon
pip install sphinxcontrib-programoutput
pip install sphinxcontrib.autorun
