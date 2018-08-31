.. _coding-guide:

Coding guide
============

.. contents::
   :backlinks: top

Committing changes
------------------

When committing changes to AbiPy, there are a few things to bear in mind.

* if your changes are non-trivial, please make an entry in the :file:`CHANGELOG.rst`
  Note that the changelog is written following the format employed by the 
  `releases <https://github.com/bitprophet/releases>`_ Sphinx extension.

* if you change the API, please document the modifications in the docstring with the ``versionadded`` role::

    .. versionadded:: 0.2
       Add new argument ``foobar``

* Are your changes python2.7 compatible?

* Can you pass the automatic tests? 

* Can you add a test to test your changes?

* If you have added new files or directories, or reorganized existing
  ones, are the new files included in the match patterns in :file:`MANIFEST.in`.  
  This file determines what goes into the source distribution of the build.

Importing and name spaces
-------------------------

For numpy_, use::

  import numpy as np
  a = np.array([1,2,3])

For matplotlib_, **avoid** using the high-level interface such as in::

  import matplotlib.pyplot as plt
  plt.plot(x, y)

and use the object-oriented API provided by |matplotlib-Axes|::

    from abipy.tools.plotting import get_ax_fig_plt
    ax, fig, plt = get_ax_fig_plt(ax=None)
    ax.plot(x, y)

It is recommended to pass the Axes ``ax`` to the plotting method and 
use the decorator ``add_fig_kwargs``::

    @add_fig_kwargs
    def plot(self, ax=None, **kwargs):
        """
        Plot the object ...

        Args:
            ax: |matplotlib-Axes| or None if a new figure should be created.

        Returns: |matplotlib-Figure|
        """
        ax, fig, plt = get_ax_fig_plt(ax=ax)
        ax.plot(self.xvals, self.yvalsm, **kwargs)
        return fig

Naming, spacing, and formatting conventions
-------------------------------------------

In general, we want to stay as closely as possible to the standard
coding guidelines for python written by Guido van Rossum in `PEP0008 <http://www.python.org/dev/peps/pep-0008>`_.

* functions and class methods: ``lower`` or ``lower_underscore_separated``
* attributes and variables: ``lower`` 
* classes: ``Upper`` or ``MixedCase``

Prefer the shortest names that are still readable.

Configure your editor to use spaces, not hard tabs. 
The standard indentation unit is always four spaces; 
if there is a file with tabs or a different number of spaces it is a bug -- please fix it.

Keep docstrings uniformly indented as in the example below, with nothing to the left of the triple quotes.  

Limit line length to (around) 90 characters. 
If you wonder why we are violating pep8 that specifies a maximum line length of 79 characters,
check out this video by Raymond Hettinger:

.. youtube:: wf-BqAjZb8M

If a logical line needs to be longer, use parentheses to break it; do not use an escaped newline.
It may be preferable to use a temporary variable to replace a single
long line with two shorter and more readable lines.

Please do not commit lines with trailing white space, as it causes noise in diffs.  

Writing examples
----------------

We have examples in subdirectories of :file:`abipy/examples`, and these are automatically
generated when the website is built to show up both in the :file:`examples`
and :file:`gallery` sections of the website.  
Many people find these examples from the website, and do not have ready access to the 
:file:`examples` directory in which they reside.  
Thus any example data that is required for the example should be added to the :file:`abipy/data` directory

Testing
-------

Abipy has a testing infrastructure based on :mod:`unittest` and pytest_.

Common test support is provided by :mod:`abipy.core.testing`, 
data files are stored in :file:`abipy/data`, in particular in :file:`abipy/data/refs` that
contains several output files that can be used for writing unit tests and examples.
