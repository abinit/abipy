#!/usr/bin/env python
"""Deployment file to facilitate releases of abipy."""
from __future__ import division, print_function

#import glob
#import os

from fabric.api import local, lcd

def test():
    """Run the unit tests."""
    local("nosetests")

def commit():
    """Commit changes."""
    local("git commit -a")

def push():
    """Push to github."""
    local("git push")

def deploy():
    """Unit tests + git commit + git push."""
    test()
    commit()
    push()

#def publish():
#    local("python setup.py release")

def makedoc():
    """Build the docs."""

    # Generate html pages from the notebooks.
    with lcd("abipy/examples/notebooks"):
       local("make html") # "ipython nbconvert --to html *.ipynb"
       #local("mv *.html ../../docs/_static")

    # Run sphynx.
    with lcd("docs"):
        local("make clean")
        local("make html")
        #local("cp _static/* _build/html/_static")

def makedoc_check():
    """Check the docstrings and the links."""
    with lcd("docs"):
        #local("make doctest")
        local("make linkcheck")

#def update_doc():
#    makedoc()
#    with lcd("docs/_build/html/"):
#        local("git add .")
#        local("git commit -a -m \"Update dev docs\"")
#        local("git push origin gh-pages")

#def release():
#    setver()
#    test()
#    publish()
#    log_ver()
#    update_doc()
