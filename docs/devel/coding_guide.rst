.. _coding-guide:

************
Coding guide
************

Committing changes
==================

When committing changes to abipy, there are a few things to bear in mind.

* if your changes are non-trivial, please make an entry in the :file:`CHANGELOG`

* if you change the API, please document it in :file:`doc/api/api_changes.rst`,

* Are your changes python2.7 compatible?

* Can you pass the automatic tests? (:file:`README-tests.rst` describes how to run the tests)

* Can you add a test to :file:`abipy/tests` to test your changes?

* if you have added new files or directories, or reorganized existing
  ones, are the new files included in the match patterns in :file:`MANIFEST.in`.  
  This file determines what goes into the source distribution of the build.


Style guide
===========

Importing and name spaces
-------------------------

For `numpy <http://www.numpy.org>`_, use::

  import numpy as np
  a = np.array([1,2,3])

For `matplotlib <http://matplotlib.org/>`_, use::

  import matplotlib as mpl
  mpl.rcParams['xtick.major.pad'] = 6

  import matplotlib.pyplot as plt
  plt.plot(x, y)

Naming, spacing, and formatting conventions
-------------------------------------------

In general, we want to stay as closely as possible to the standard
coding guidelines for python written by Guido in 
`PEP0008 <http://www.python.org/dev/peps/pep-0008>`_.

* functions and class methods: ``lower`` or ``lower_underscore_separated``

* attributes and variables: ``lower`` 

* classes: ``Upper`` or ``MixedCase``

Prefer the shortest names that are still readable.

Configure your editor to use spaces, not hard tabs. 
The standard indentation unit is always four spaces; 
if there is a file with tabs or a different number of spaces it is a bug -- please fix it.
To detect and fix these and other whitespace errors (see below),
use :file:`tools/reindent.py` as a command-line script.  
Unless you are sure your editor always does the right thing, please use reindent.py before committing your changes.

Keep docstrings uniformly indented as in the example below, with nothing to the left of the triple quotes.  

Limit line length to 90 characters.  
If a logical line needs to be longer, use parentheses to break it; do not use an escaped newline.
It may be preferable to use a temporary variable to replace a single
long line with two shorter and more readable lines.

Please do not commit lines with trailing white space, as it causes noise in diffs.  

Writing examples
================

We have examples in subdirectories of :file:`abipy/examples`, and these are automatically
generated when the website is built to show up both in the :file:`examples`
and :file:`gallery` sections of the website.  
Many people find these examples from the website, and do not have ready access to the 
:file:`examples` directory in which they reside.  
Thus any example data that is required for the example should be added to the :file:`example/data` directory

Testing
=======

Abipy has a testing infrastructure based on :mod:`unittest`.
The tests are in :mod:`abipy.tests`, data files are strore in :file:`abipy/tests/data`.
