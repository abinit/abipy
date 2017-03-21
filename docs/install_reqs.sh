#!/bin/bash
set -ev  # exit on first error, print each command

conda config --set always_yes yes --set changeps1 no
conda install sphinx matplotlib ipython
pip install sphinxcontrib.napoleon
pip install sphinxcontrib-programoutput
#conda install sphinx_rtd_theme
#pip install sphinxcontrib.autorun
